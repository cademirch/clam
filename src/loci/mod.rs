pub mod d4_bgzf;
pub mod d4_tasks;

use anyhow::{bail, Result};
use d4::{
    ptab::PTablePartitionWriter, stab::SecondaryTablePartWriter, Chrom, D4FileBuilder,
    D4FileMerger, D4FileWriter, Dictionary,
};
use log::{debug, warn};
use serde::Deserialize;
use std::collections::HashMap;
use std::collections::{hash_map::Entry, HashSet};
use std::fs::File;
use std::io::Write;
use std::path::{Path, PathBuf};
use tempfile::NamedTempFile;

use clap::{Parser, Subcommand, ArgGroup};
use camino::Utf8PathBuf;


#[derive(Parser, Debug, Clone)]
#[command(author, version, about = "Calculate callable sites from depth statistics.")]

pub struct LociArgs {
    /// Path to input D4 file
    #[arg(required(true))]
    pub infile: Utf8PathBuf,
    /// Path to output file. Extensions allowed: {".bed", ".d4"}. ".bed" is mutually exclusive with --populations
    #[arg(required(true))]
    pub outfile: Utf8PathBuf,
    /// Minimum depth to consider site callable per individual
    #[arg(short = 'm', long = "min-depth", default_value_t = 0.0, conflicts_with("threshold_file"))]
    pub min_depth: f64,
    /// Maximum depth to consider site callable per individual
    #[arg(short = 'M', long = "max-depth", default_value_t = f64::INFINITY, conflicts_with("threshold_file"))]
    pub max_depth: f64,
    /// Proportion of samples passing thresholds at site to consider callable
    #[arg(short = 'd', long = "depth-proportion", default_value_t = 1.0, conflicts_with("population_file"))]
    pub depth_proportion: f64,
    /// Minimum mean depth across all samples at site to consider callable
    #[arg(short = 'u', long = "min-mean-depth", default_value_t = 0.0, conflicts_with("population_file"))]
    pub mean_depth_min: f64,
    /// Maximum mean depth across all samples at site to consider callable
    #[arg(short = 'U', long = "max-mean-depth", default_value_t = f64::INFINITY, conflicts_with("population_file"))]
    pub mean_depth_max: f64,
    /// Output number of individuals callable at site.
    #[arg(short = 'c', long = "output-counts", default_value_t = false)]
    pub output_counts: bool,
    /// Number of threads to use
    #[arg(short = 't', long = "threads", default_value_t = 1)]
    pub threads: usize,
    /// Path to file that defines populations. Tab separated: sample, population_name
    #[arg(short = 'p', long = "population-file")]
    pub population_file: Option<Utf8PathBuf>,
    /// Path to file that defines per-chromosome individual level thresholds. Tab separated: chrom, min, max
    #[arg(long = "thresholds-file")]
    pub threshold_file: Option<Utf8PathBuf>,
    /// Comma separated list of chromosomes to exclude
    #[arg(short = 'x', value_delimiter = ',', num_args = 1.., conflicts_with("exclude_file"))]
    pub exclude: Option<Vec<String>>,
    /// Path to file with chromosomes to exclude, one per line
    #[arg(long = "exclude-file", conflicts_with("exclude"))]
    pub exclude_file: Option<Utf8PathBuf>,
}

pub enum D4Reader {
    D4(d4::D4TrackReader),
    Bgzf(d4_bgzf::BgzfD4MatrixReader),
}

impl D4Reader {
    pub fn chrom_regions(
        &self,
        thresholds: Thresholds,
        exclude_chrs: Option<&HashSet<String>>,
    ) -> Result<Vec<ChromRegion>> {
        match self {
            D4Reader::D4(reader) => {
                prepare_chrom_regions(reader.chrom_regions(), thresholds, exclude_chrs)
            }
            D4Reader::Bgzf(reader) => {
                prepare_chrom_regions(reader.chrom_regions(), thresholds, exclude_chrs)
            }
        }
    }

    pub fn run_tasks(
        &self,
        loci_args: &LociArgs,
        chrom_regions: Vec<ChromRegion>,
        sample_refs: Option<Vec<String>>, // Updated to Vec<String>
    ) -> Result<Vec<(String, u32, Vec<CallableRegion>)>> {
        match self {
            D4Reader::Bgzf(_) => d4_bgzf::run_bgzf_tasks(
                loci_args.infile.clone(),
                loci_args.threads,
                sample_refs.clone(), // Pass Vec<String> directly
                chrom_regions,
                (loci_args.mean_depth_min, loci_args.mean_depth_max),
                loci_args.depth_proportion,
                loci_args.output_counts,
            ),
            D4Reader::D4(_) => {
                let tracks = d4_tasks::prepare_tracks_from_file(loci_args.infile.clone(), sample_refs)?;
                d4_tasks::run_tasks_on_tracks(
                    tracks,
                    chrom_regions,
                    (loci_args.mean_depth_min, loci_args.mean_depth_max),
                    loci_args.depth_proportion,
                    loci_args.output_counts,
                )
            }
        }
    }
    
}

#[derive(Clone)]
pub enum Thresholds {
    Fixed((f64, f64)),                          // For single threshold tuple
    PerChromosome(HashMap<String, (f64, f64)>), // For chromosome-specific thresholds
}

#[derive(Debug, Deserialize)]
struct ThresholdRecord {
    chrom: String,
    min: f64,
    max: f64,
}

pub fn read_threshold_file<P: AsRef<Path>>(
    threshold_file: P,
) -> Result<HashMap<String, (f64, f64)>> {
    let file = File::open(&threshold_file).expect(&format!(
        "Failed to open population file: {}",
        threshold_file.as_ref().display()
    ));
    let mut reader = csv::ReaderBuilder::new().delimiter(b'\t').from_reader(file);
    let mut filter_map: HashMap<String, (f64, f64)> = HashMap::new();

    for result in reader.deserialize() {
        let record: ThresholdRecord = result?;

        match filter_map.entry(record.chrom) {
            Entry::Vacant(entry) => {
                entry.insert((record.min, record.max));
            }
            Entry::Occupied(mut entry) => {
                warn!(
                    "Duplicate chromosome in thresholds file! {} already exists with min: {}, max: {}. Overwriting with new values min: {}, max: {}.",
                    entry.key(),
                    entry.get().0,
                    entry.get().1,
                    record.min,
                    record.max
                );
                entry.insert((record.min, record.max));
            }
        }
    }

    Ok(filter_map)
}

#[derive(Clone, Debug, PartialEq)]
pub struct CallableRegion {
    pub count: u32,
    pub begin: u32,
    pub end: u32,
}

pub fn merge_d4_files<P: AsRef<Path>>(outpath: P, input: Vec<P>, names: Vec<&str>) -> Result<()> {
    let mut merger = D4FileMerger::new(outpath);
    for (file, name) in input.iter().zip(names.iter()) {
        merger = merger.add_input_with_tag(file, name)
    }
    merger.merge()?;
    Ok(())
}

pub fn write_d4<P: AsRef<Path>>(
    regions: Vec<(String, u32, Vec<CallableRegion>)>, // Generic input structure
    chroms: Vec<Chrom>,
    output_path: Option<P>,
) -> Result<PathBuf> {
    // Initialize the D4 file writer

    let output_path: PathBuf = if let Some(output_path) = output_path {
        output_path.as_ref().to_path_buf() // Convert P to PathBuf
    } else {
        let temp_file = NamedTempFile::new()?;
        temp_file.into_temp_path().to_path_buf() // Convert NamedTempFile path to PathBuf
    };

    let mut builder = D4FileBuilder::new(output_path.to_path_buf());
    builder.append_chrom(chroms.into_iter());
    builder.set_denominator(1.0); // Set denominator (you can adjust as needed)
    builder.set_dictionary(Dictionary::new_simple_range_dict(0, 1)?); // Set dictionary for encoding

    // Create the D4 file writer
    let mut d4_writer: D4FileWriter = builder.create()?;

    let mut partitions = d4_writer.parallel_parts(None)?;
    for (idx, partition) in partitions.iter().enumerate() {
        let region_info = partition.0.region();
        debug!(
            "Partition {}: chrom = {}, region = {}-{}",
            idx, region_info.0, region_info.1, region_info.2
        );
    }

    // Process each tuple in `regions`, which contains a chromosome, starting position, and callable regions
    for (chrom, begin, callable_regions) in regions {
        let mut current_partition = 0; // Track the current partition

        // Process each CallableRegion
        for region in callable_regions.into_iter() {
            let start = region.begin + begin;
            let end = region.end + begin;

            let mut pos = start;
            while pos < end {
                let region_info = partitions[current_partition].0.region();

                // Check if we are in the correct region for this partition
                if region_info.0 != chrom || region_info.1 < pos || region_info.2 <= pos {
                    // If not, find the correct partition
                    if let Some((idx, _)) = partitions.iter().enumerate().find(|(_, part)| {
                        let reg = part.0.region();
                        reg.0 == chrom && reg.1 <= pos && pos <= reg.2
                    }) {
                        current_partition = idx;
                    } else {
                        debug!(
                            "Could not find matching partition for chrom: {}, pos: {}",
                            chrom, pos
                        );
                        continue;
                    }
                }

                let record_end = region_info.2.min(end);

                let mut primary_encoder = partitions[current_partition].0.make_encoder();
                let mut secondary_table = &mut partitions[current_partition].1;

                for current_pos in pos..record_end {
                    if !primary_encoder.encode(current_pos as usize, region.count as i32) {
                        debug!("Encoding secondary table for pos: {}", current_pos); // Debugging output
                        secondary_table.encode(current_pos as u32, region.count as i32)?;
                    }
                }

                pos = record_end; // Move the position to the next unencoded part
            }

            debug!(
                "Flushing secondary table for partition {}",
                current_partition
            ); // Debugging output
               // Flush the current partition after processing
            partitions[current_partition].1.flush()?;
        }
    }

    // Finish writing all partitions
    for (idx, (_, mut secondary_table)) in partitions.into_iter().enumerate() {
        debug!("Finishing partition {}", idx); // Debugging output
        secondary_table.finish()?;
    }

    debug!("D4 file writing complete."); // Final debug output
    Ok(output_path)
}

pub fn write_bed<P: AsRef<Path>>(
    output_path: P,
    regions: Vec<(String, u32, Vec<CallableRegion>)>,
    output_counts: bool,
) -> Result<()> {
    let mut file = File::create(output_path)?;

    for (chrom, begin, callable_regions) in regions {
        for region in callable_regions.into_iter() {
            if output_counts {
                writeln!(
                    file,
                    "{}\t{}\t{}\t{}",
                    chrom,
                    region.begin + begin,
                    region.end + begin,
                    region.count
                )?;
            } else {
                writeln!(
                    file,
                    "{}\t{}\t{}",
                    chrom,
                    region.begin + begin,
                    region.end + begin
                )?;
            }
        }
    }

    Ok(())
}

#[derive(Clone)]
pub struct ChromRegion {
    pub chr: String,
    pub begin: u32,
    pub end: u32,
    pub min_filter: f64,
    pub max_filter: f64,
}

pub fn prepare_chrom_regions(
    chrom_regions: Vec<(&str, u32, u32)>,
    thresholds: Thresholds,
    exclude_chrs: Option<&HashSet<String>>,
) -> Result<Vec<ChromRegion>> {
    // Apply thresholds and exclusion logic
    let chrom_filters: Vec<ChromRegion> = match thresholds {
        Thresholds::Fixed(thresh) => chrom_regions
            .into_iter()
            .map(|(chr, start, end)| ChromRegion {
                chr: chr.to_string(),
                begin: start,
                end,
                min_filter: thresh.0,
                max_filter: thresh.1,
            })
            .collect(),

        Thresholds::PerChromosome(ref filter_map) => {
            let chroms_with_filters = add_filters_to_chroms(
                chrom_regions
                    .iter()
                    .map(|(chr, start, end)| (*chr, *start, *end))
                    .collect(),
                filter_map.clone(),
            )?;
            chroms_with_filters
                .into_iter()
                .map(|(chr, start, end, min_filter, max_filter)| ChromRegion {
                    chr,
                    begin: start,
                    end,
                    min_filter,
                    max_filter,
                })
                .collect()
        }
    };

    Ok(chrom_filters
        .into_iter()
        .filter(|region| {
            if let Some(exclude_set) = exclude_chrs {
                !exclude_set.contains(&region.chr)
            } else {
                true
            }
        })
        .collect())
}

pub fn add_filters_to_chroms(
    d4_chrom_regions: Vec<(&str, u32, u32)>,
    filter_map: HashMap<String, (f64, f64)>,
) -> Result<Vec<(String, u32, u32, f64, f64)>> {
    let mut result = vec![];
    for (chrom, start, end) in &d4_chrom_regions {
        if let Some(&(min_filter, max_filter)) = filter_map.get(*chrom) {
            result.push((chrom.to_string(), *start, *end, min_filter, max_filter));
        } else {
            // Issue a warning if no filter is found and apply default filters
            warn!(
                "No filter found for chromosome '{}', using default filters (0.0, f64::INFINITY).",
                chrom
            );
            result.push((chrom.to_string(), *start, *end, 0.0, f64::INFINITY));
        }
    }

    // Ensure every filter in chrom_filters was used; if not, bail
    for (chrom, (_, _)) in filter_map.iter() {
        if !d4_chrom_regions
            .iter()
            .any(|(region_chrom, _, _)| region_chrom == chrom)
        {
            bail!(
                "Filter provided for chromosome '{}' not found in D4 chrom regions.",
                chrom
            );
        }
    }

    Ok(result)
}

#[cfg(test)]
mod tests {

    use super::*;
    use log::debug;
    use std::collections::{HashMap, HashSet};
    fn init() {
        let _ = env_logger::builder()
            .target(env_logger::Target::Stdout)
            .filter_level(log::LevelFilter::Trace)
            .is_test(true)
            .try_init();
    }

    #[test]
    fn test_run_and_write_d4() -> Result<()> {
        // Initialize logger for debugging in tests
        init();
        debug!("Starting run and write test...");
        // Step 1: Load the D4 file and read tracks
        let d4_file_path = "tests/data/merged.d4"; // Path to your D4 file
        let samples: Vec<String> = vec!["0.per".to_string(), "2.per".to_string()];

        let tracks = d4_tasks::prepare_tracks_from_file(d4_file_path, Some(samples.clone()))?;
        // number of tracks should be same as number of samples specified
        assert_eq!(samples.len(), tracks.len());
        let thresholds = Thresholds::Fixed((0.0, f64::INFINITY)); // Default thresholds
        let mean_thresholds = (0.0, f64::INFINITY); // Example mean thresholds
        let depth_proportion = 0.0; // Default depth proportion
        let output_counts = true; // Default setting, adjust if needed
        let exclude_chrs: Option<HashSet<String>> = None; // No chromosomes to exclude by default

        let rdr: d4::D4TrackReader = d4::D4TrackReader::open_first_track(d4_file_path)?;
        let chrom_regions =
            super::prepare_chrom_regions(rdr.chrom_regions(), thresholds, exclude_chrs.as_ref())?;
        let chroms: Vec<d4::Chrom> = chrom_regions
            .clone()
            .into_iter()
            .map(|r| d4::Chrom {
                name: r.chr.to_string(),
                size: r.end.try_into().unwrap(),
            })
            .collect();
        drop(rdr);

        let task_output = d4_tasks::run_tasks_on_tracks(
            tracks,
            chrom_regions,
            mean_thresholds,
            depth_proportion,
            output_counts,
        )?;

        debug!("Writing d4 file...");
        let output_path = write_d4::<PathBuf>(task_output, chroms, None)?;

        assert!(
            std::fs::metadata(output_path.to_path_buf()).is_ok(),
            "Output D4 file should be created"
        );

        std::fs::remove_file(output_path).expect("Failed to clean up output file");

        debug!("Successfully processed and wrote output to D4.");
        Ok(())
    }

    #[test]
    fn test_add_filters_to_chroms_success() {
        // Setup input
        let d4_chrom_regions = vec![
            ("chr1", 0, 1000),
            ("chr2", 1000, 2000),
            ("chr3", 2000, 3000),
        ];
        let mut filter_map = HashMap::new();
        filter_map.insert("chr1".to_string(), (0.5, 1.5));
        filter_map.insert("chr2".to_string(), (0.1, 0.9));
        filter_map.insert("chr3".to_string(), (0.2, 1.0));

        // Expected output
        let expected = vec![
            ("chr1".to_string(), 0, 1000, 0.5, 1.5),
            ("chr2".to_string(), 1000, 2000, 0.1, 0.9),
            ("chr3".to_string(), 2000, 3000, 0.2, 1.0),
        ];

        // Call function
        let result = add_filters_to_chroms(d4_chrom_regions, filter_map).unwrap();

        // Assert the result matches expected
        assert_eq!(result, expected);
    }

    #[test]
    fn test_add_filters_to_chroms_missing_filter() {
        // Setup input with one chromosome missing from the filter map
        let d4_chrom_regions = vec![
            ("chr1", 0, 1000),
            ("chr2", 1000, 2000),
            ("chr3", 2000, 3000),
        ];
        let mut filter_map = HashMap::new();
        filter_map.insert("chr1".to_string(), (0.5, 1.5));
        filter_map.insert("chr2".to_string(), (0.1, 0.9));
        // chr3 is missing in the filter map

        // Expected output
        let expected = vec![
            ("chr1".to_string(), 0, 1000, 0.5, 1.5),
            ("chr2".to_string(), 1000, 2000, 0.1, 0.9),
            ("chr3".to_string(), 2000, 3000, 0.0, f64::INFINITY), // Default filter for chr3
        ];

        // Call function
        let result = add_filters_to_chroms(d4_chrom_regions, filter_map).unwrap();

        // Assert the result matches expected
        assert_eq!(result, expected);
    }

    #[test]
    fn test_add_filters_to_chroms_unmatched_filter() {
        // Setup input with an extra chromosome in the filter map
        let d4_chrom_regions = vec![("chr1", 0, 1000), ("chr2", 1000, 2000)];
        let mut filter_map = HashMap::new();
        filter_map.insert("chr1".to_string(), (0.5, 1.5));
        filter_map.insert("chr2".to_string(), (0.1, 0.9));
        filter_map.insert("chrX".to_string(), (0.3, 1.2)); // Unmatched chromosome filter

        // Call function, expecting an error
        let result = add_filters_to_chroms(d4_chrom_regions, filter_map);

        // Assert the function returns an error due to unmatched chromosome in the filter map
        assert!(result.is_err());
        let err_message = result.unwrap_err().to_string();
        assert!(err_message
            .contains("Filter provided for chromosome 'chrX' not found in D4 chrom regions."));
    }
}
