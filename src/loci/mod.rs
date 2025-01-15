pub mod d4_bgzf;
pub mod d4_tasks;

use std::collections::hash_map::Entry;
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::Write;
use std::path::{Path, PathBuf};

use anyhow::{bail, Result, Context};
use camino::Utf8PathBuf;
use clap::Parser;
use d4::index::D4IndexCollection;
use d4::ptab::PTablePartitionWriter;
use d4::stab::SecondaryTablePartWriter;
use d4::{Chrom, D4FileBuilder, D4FileMerger, D4FileWriter, Dictionary};
use log::{debug, trace, warn};
use rayon::prelude::*;
use serde::Deserialize;
use tempfile::NamedTempFile;

#[derive(Parser, Debug, Clone)]
#[command(
    author,
    version,
    about = "Calculate callable sites from depth statistics."
)]

pub struct LociArgs {
    /// Path to input D4 file
    #[arg(required(true))]
    pub infile: Utf8PathBuf,

    /// Output file prefix. The extension will be added automatically based on the `--no-counts` flag.
    #[arg(required(true))]
    pub outprefix: Utf8PathBuf,

    /// Minimum depth to consider site callable per individual
    #[arg(
        short = 'm',
        long = "min-depth",
        default_value_t = 0.0,
        conflicts_with("threshold_file")
    )]
    pub min_depth: f64,

    /// Maximum depth to consider site callable per individual
    #[arg(short = 'M', long = "max-depth", default_value_t = f64::INFINITY, conflicts_with("threshold_file"))]
    pub max_depth: f64,

    /// Proportion of samples passing thresholds at site to consider callable. Ignored when outputting counts.
    #[arg(
        short = 'd',
        long = "depth-proportion",
        default_value_t = 0.0,
        conflicts_with("population_file")
    )]
    pub depth_proportion: f64,

    /// Minimum mean depth across all samples at site to consider callable. Ignored when outputting counts.
    #[arg(
        short = 'u',
        long = "min-mean-depth",
        default_value_t = 0.0,
        conflicts_with("population_file")
    )]
    pub mean_depth_min: f64,

    /// Maximum mean depth across all samples at site to consider callable. Ignored when outputting counts.
    #[arg(short = 'U', long = "max-mean-depth", default_value_t = f64::INFINITY, conflicts_with("population_file"))]
    pub mean_depth_max: f64,

    /// Disable outputting counts; produces a .bed file instead.
    #[arg(
        long = "no-counts",
        default_value_t = false,
        conflicts_with("population_file")
    )]
    pub no_counts: bool,

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

impl LociArgs {
    /// Resolve the final output file path with the appropriate extension based on `--no-counts` flag.
    pub fn resolve_output_file(&self) -> Utf8PathBuf {
        let mut output_file = self.outprefix.clone();
        let extension = if self.no_counts { "bed" } else { "d4" };
        output_file.set_extension(extension);
        output_file
    }
}

pub enum D4Reader {
    D4(d4::D4TrackReader),
    Bgzf(d4_bgzf::BGZID4MatrixReader),
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
        sample_refs: Option<Vec<String>>, 
    ) -> Result<Vec<(String, u32, Vec<CallableRegion>)>> {
        let output_counts = !loci_args.no_counts;
        debug!("output_counts: {output_counts}");
        match self {
            D4Reader::Bgzf(_) => d4_bgzf::run_bgzf_tasks(
                loci_args.infile.clone(),
                loci_args.threads,
                sample_refs.clone(), // Pass Vec<String> directly
                chrom_regions,
                (loci_args.mean_depth_min, loci_args.mean_depth_max),
                loci_args.depth_proportion,
                output_counts,
            ),
            D4Reader::D4(_) => {
                let tracks =
                    d4_tasks::prepare_tracks_from_file(loci_args.infile.clone(), sample_refs)?;
                d4_tasks::run_tasks_on_tracks(
                    tracks,
                    chrom_regions,
                    (loci_args.mean_depth_min, loci_args.mean_depth_max),
                    loci_args.depth_proportion,
                    output_counts,
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
        "Failed to open threshold file: {}",
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

pub fn write_d4_parallel<P: AsRef<Path>>(
    regions: Vec<(String, u32, Vec<CallableRegion>)>,
    chroms: Vec<Chrom>,
    output_path: Option<P>,
) -> Result<PathBuf> {
    // Initialize the D4 file writer
    let output_path: PathBuf = if let Some(output_path) = output_path {
        output_path.as_ref().to_path_buf()
    } else {
        let temp_file = NamedTempFile::new()?;
        temp_file.into_temp_path().to_path_buf()
    };

    let mut builder = D4FileBuilder::new(output_path.clone());
    builder.append_chrom(chroms.into_iter());
    builder.set_denominator(1.0);
    builder.set_dictionary(Dictionary::new_simple_range_dict(0, 1)?);

    // Create the D4 file writer
    let mut d4_writer: D4FileWriter = builder.create()?;
    let partitions = d4_writer.parallel_parts(Some(10_000_000))?;
    trace!("Created {} partitions for writing", partitions.len());

    let partitions_result: Result<()> = partitions
        .into_par_iter()
        .map(|(mut partition, mut secondary_table)| {
            let (chrom, left, right) = partition.region();
            let chrom = chrom.to_string();
            let mut primary_encoder = partition.make_encoder();
            let mut last = left;

            // Encapsulate encoding logic in a closure
            let mut write_value = |pos: u32, value: i32| {
                if !primary_encoder.encode(pos as usize, value) {
                    secondary_table.encode(pos, value).unwrap();
                    // trace!("Secondary encode: Wrote {} at {}:{}", value, chrom, pos);
                } else {
                    // trace!("Primary encode: Wrote {} at {}:{}", value, chrom, pos);
                }
            };

            // Process regions within the current partition
            if let Some((_, begin_offset, region_list)) = regions
                .iter()
                .find(|(region_chrom, _, _)| *region_chrom == chrom)
            {
                for region in region_list {
                    let start = region.begin + *begin_offset;
                    let end = region.end + *begin_offset;

                    if end <= left || start >= right {
                        // Skip regions outside of the partition bounds
                        continue;
                    }

                    // Clip the region to fit within the partition bounds
                    let clipped_start = start.max(left);
                    let clipped_end = end.min(right);

                    for pos in last..clipped_start {
                        write_value(pos, 0);
                    }

                    for pos in clipped_start..clipped_end {
                        write_value(pos, region.count as i32);
                    }

                    last = clipped_end;
                }
            }

            // Fill remaining positions with zero
            for pos in last..right {
                write_value(pos, 0);
            }

            // Flush the secondary table
            secondary_table.flush()?;
            secondary_table.finish()?;
            Ok(())
        })
        .collect();

    // Check for any errors during the parallel operation
    if let Err(e) = partitions_result {
        return Err(e.into());
    }
    
    drop(d4_writer);
    log::info!("D4 file writing complete.");
    
    let mut index_collection = D4IndexCollection::open_for_write(output_path.clone())
        .with_context(|| {
            format!(
                "Failed to open D4 index collection for writing at {:?}",
                output_path
            )
        })?;

    index_collection
        .create_secondary_frame_index()
        .with_context(|| "Failed to create secondary frame index for the D4 index collection")?;
    Ok(output_path)
}

pub fn write_bed<P: AsRef<Path>>(
    output_path: P,
    regions: Vec<(String, u32, Vec<CallableRegion>)>,
) -> Result<()> {
    let mut file = File::create(output_path)?;

    for (chrom, begin, callable_regions) in regions {
        for region in callable_regions.into_iter() {
            writeln!(
                file,
                "{}\t{}\t{}",
                chrom,
                region.begin + begin,
                region.end + begin
            )?;
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

    use std::collections::HashMap;

    use super::*;
    fn init() {
        let _ = env_logger::builder()
            .target(env_logger::Target::Stdout)
            .filter_level(log::LevelFilter::Trace)
            .is_test(true)
            .try_init();
    }

    #[test]
    fn test_add_filters_to_chroms_success() {
        // Setup input
        init();
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
