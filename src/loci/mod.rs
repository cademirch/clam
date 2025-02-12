pub mod counts;
pub mod d4_bgzf;
pub mod d4_tasks;
pub mod gvcf;
pub mod io;
pub mod regions;
pub mod thresholds;

use std::collections::hash_map::Entry;
use std::collections::{HashMap, HashSet, VecDeque};
use std::fs::File;
use std::io::{BufRead, BufReader, Write};
use std::path::{Path, PathBuf};

use anyhow::{bail, Context, Result};
use camino::Utf8PathBuf;
use clap::Parser;
use d4::index::D4IndexCollection;
use d4::ptab::PTablePartitionWriter;
use d4::stab::SecondaryTablePartWriter;
use d4::{Chrom, D4FileBuilder, D4FileMerger, D4FileWriter, Dictionary};
use indicatif::ProgressBar;
use log::{debug, trace, warn};
use rayon::prelude::*;
use serde::Deserialize;
use tempfile::NamedTempFile;

#[derive(Debug, Clone, PartialEq)]
pub enum InputMode {
    MergedD4(PathBuf),            // Single merged D4 file
    MultiD4(VecDeque<PathBuf>),   // Multiple per-sample D4 files
    MultiGVCF(VecDeque<PathBuf>), // Multiple per-sample GVCF files
}

// struct ProcessConfig {
//     per_sample_thresholds: thresholds::Thresholds,
//     mean_thresholds: thresholds::Thresholds,
//     exclude_chrs: Option<HashSet<String>>,
//     outdir: PathBuf,
//     depth_prop: f64,
//     population_map: Option<super::utils::PopulationMapping>,
// }

#[derive(Parser, Debug, Clone)]
#[command(
    author,
    version,
    about = "Calculate callable sites from depth statistics."
)]
pub struct LociArgs {
    /// Input files (D4 format by default)
    #[arg(required_unless_present = "filelist")]
    pub input: Vec<Utf8PathBuf>,

    /// Path to file containing list of input files, one per line
    #[arg(short = 'f', long = "filelist")]
    pub filelist: Option<Utf8PathBuf>,

    /// Use GVCF format instead of D4
    #[arg(long = "gvcf")]
    pub gvcf: bool,

    /// Input is a merged D4 file
    #[arg(long = "merged", conflicts_with = "gvcf")]
    pub merged: bool,

    /// Output directory
    #[arg(short = 'o', required(true))]
    pub outdir: Utf8PathBuf,

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

    /// Proportion of samples passing thresholds at site to consider callable
    #[arg(
        short = 'd',
        long = "depth-proportion",
        default_value_t = 0.0,
        conflicts_with("population_file")
    )]
    pub depth_proportion: f64,

    /// Minimum mean depth across all samples at site to consider callable
    #[arg(
        short = 'u',
        long = "min-mean-depth",
        default_value_t = 0.0,
        conflicts_with("population_file")
    )]
    pub mean_depth_min: f64,

    /// Maximum mean depth across all samples at site to consider callable
    #[arg(short = 'U', long = "max-mean-depth", default_value_t = f64::INFINITY, conflicts_with("population_file"))]
    pub mean_depth_max: f64,

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

    // Fields that get initialized
    #[arg(skip)]
    pub exclude_chrs: Option<HashSet<String>>,

    #[arg(skip)]
    pub population_map: Option<super::utils::PopulationMapping>,
}

impl LociArgs {
    /// Initialize processing-related fields and set up environment
    pub fn initialize(&mut self) -> Result<()> {
        // Set up thread pool
        rayon::ThreadPoolBuilder::new()
            .num_threads(self.threads)
            .build_global()?;

        // Setup output directory
        if !self.outdir.exists() {
            std::fs::create_dir_all(&self.outdir)?;
        }

        // Initialize excluded chromosomes
        self.exclude_chrs =
            super::utils::get_exclude_chromosomes(&self.exclude, &self.exclude_file)?;

        // Initialize population mapping if provided
        if let Some(pop_file) = &self.population_file {
            let file = File::open(pop_file).with_context(|| {
                format!(
                    "Failed to open population file at path: {}",
                    pop_file.as_str()
                )
            })?;
            self.population_map = Some(super::utils::PopulationMapping::from_path(file, None)?);
        }

        Ok(())
    }

    pub fn get_per_sample_thresholds(&self) -> thresholds::Thresholds {
        self.threshold_file.as_ref().map_or_else(
            || thresholds::Thresholds::Fixed((self.min_depth, self.max_depth)),
            |file| {
                thresholds::read_threshold_file(file)
                    .map(thresholds::Thresholds::PerChromosome)
                    .unwrap_or_else(|_| {
                        thresholds::Thresholds::Fixed((self.min_depth, self.max_depth))
                    })
            },
        )
    }

    pub fn get_mean_thresholds(&self) -> thresholds::Thresholds {
        thresholds::Thresholds::Fixed((self.mean_depth_min, self.mean_depth_max))
    }

    pub fn determine_input_mode(&self) -> Result<InputMode> {
        let input_files = self.get_input_files(None)?;

        if input_files.is_empty() {
            anyhow::bail!("At least one input file is required");
        }

        match (self.gvcf, self.merged) {
            (true, _) => Ok(InputMode::MultiGVCF(input_files)),
            (false, true) => {
                if input_files.len() > 1 {
                    anyhow::bail!("Merged D4 mode expects exactly one input file");
                }
                Ok(InputMode::MergedD4(input_files.into_iter().next().unwrap()))
            }
            (false, false) => Ok(InputMode::MultiD4(input_files)),
        }
    }

    /// Get all input files from either direct input or filelist
    pub fn get_input_files(&self, patterns: Option<&[String]>) -> Result<VecDeque<PathBuf>> {
        let mut paths = VecDeque::new();

        // Handle filelist if provided
        if let Some(filelist) = &self.filelist {
            let file = File::open(filelist)
                .with_context(|| format!("Failed to open filelist: {}", filelist))?;

            let reader = BufReader::new(file);

            for line in reader.lines() {
                let path = line
                    .with_context(|| format!("Failed to read line from filelist: {}", filelist))?;

                // Skip empty lines and comments
                if path.trim().is_empty() || path.starts_with('#') {
                    continue;
                }

                let path = PathBuf::from(path);

                // Check if path exists
                if !path.exists() {
                    anyhow::bail!("File does not exist: {}", path.display());
                }

                // If patterns are provided, check if any match
                if let Some(patterns) = patterns {
                    let path_str = path.to_string_lossy();
                    if !patterns.iter().any(|pattern| path_str.contains(pattern)) {
                        continue;
                    }
                }

                paths.push_back(path);
            }
        } else {
            // Handle direct input files
            for path in &self.input {
                let path = PathBuf::from(path);

                // Check if path exists
                if !path.exists() {
                    anyhow::bail!("File does not exist: {}", path.display());
                }

                // If patterns are provided, check if any match
                if let Some(patterns) = patterns {
                    let path_str = path.to_string_lossy();
                    if !patterns.iter().any(|pattern| path_str.contains(pattern)) {
                        continue;
                    }
                }

                paths.push_back(path);
            }
        }

        if paths.is_empty() {
            anyhow::bail!("No valid input files found");
        }

        Ok(paths)
    }
}

pub fn process(mut args: LociArgs, progress_bar: Option<ProgressBar>) -> Result<()> {
    args.initialize()?;

    match args.determine_input_mode()? {
        InputMode::MergedD4(file) => process_merged_d4(file, &args, progress_bar),
        InputMode::MultiD4(files) => process_multi_d4(files, &args, progress_bar),
        InputMode::MultiGVCF(files) => process_gvcf(files, &args, progress_bar),
    }
}

fn process_gvcf(
    files: VecDeque<PathBuf>,
    args: &LociArgs,
    progress_bar: Option<ProgressBar>,
) -> Result<()> {
    let chroms: Vec<d4::Chrom> = {
        let reader = gvcf::GvcfReader::from_path(files.front().unwrap())?;
        let chrom_regions = reader.chrom_regions();
        chrom_regions
            .iter()
            .map(|r| d4::Chrom {
                name: r.0.to_string(),
                size: r.2.try_into().unwrap(),
            })
            .collect()
    };

    if let Some(population_map) = &args.population_map {
        let progress_bar = if let Some(bar) = progress_bar {
            bar.set_length(files.len() as u64);
            Some(bar)
        } else {
            None
        };
        let mut temp_file_paths = Vec::with_capacity(population_map.num_populations());
        for (idx, samples) in population_map
            .get_samples_per_population()
            .iter()
            .enumerate()
        {
            let sample_names: Vec<String> = samples.iter().map(|&s| s.to_string()).collect();

            let res = d4_bgzf::run_tasks(
                files.clone(),
                Some(&sample_names),
                args.clone(),
                progress_bar.clone(),
            )?;

            temp_file_paths.push(io::write_d4_parallel::<PathBuf>(
                &res,
                chroms.clone(),
                None,
            )?);

            let population_name = population_map.get_popname_refs()[idx];
            io::write_bed(
                args.outdir
                    .join(format!("{}_callable_sites.bed", population_name))
                    .as_std_path()
                    .to_path_buf(),
                &res,
            );
        }

        io::merge_d4_files(
            args.outdir
                .join("callable_sites.d4")
                .as_std_path()
                .to_path_buf(),
            temp_file_paths,
            population_map.get_popname_refs(),
        )?;
    } else {
        let res = d4_bgzf::run_tasks(files, None, args.clone(), progress_bar)?;
        io::write_d4_parallel::<PathBuf>(
            &res,
            chroms,
            Some(
                args.outdir
                    .join("callable_sites.d4")
                    .as_std_path()
                    .to_path_buf(),
            ),
        )?;
        io::write_bed(
            args.outdir
                .join("callable_sites.bed")
                .as_std_path()
                .to_path_buf(),
            &res,
        )?;
    }

    Ok(())
}

fn process_multi_d4(
    files: VecDeque<PathBuf>,
    args: &LociArgs,
    progress_bar: Option<ProgressBar>,
) -> Result<()> {
    let reader = d4_bgzf::BGZID4TrackReader::from_path(files.front().unwrap(), None)?;
    let chrom_regions = reader.chrom_regions();
    let chroms: Vec<d4::Chrom> = chrom_regions
        .iter()
        .map(|r| d4::Chrom {
            name: r.0.to_string(),
            size: r.2.try_into().unwrap(),
        })
        .collect();

    if let Some(population_map) = &args.population_map {
        let progress_bar = if let Some(bar) = progress_bar {
            bar.set_length(files.len() as u64);
            Some(bar)
        } else {
            None
        };
        let mut temp_file_paths = Vec::with_capacity(population_map.num_populations());
        for (idx, samples) in population_map
            .get_samples_per_population()
            .iter()
            .enumerate()
        {
            let sample_names: Vec<String> = samples.iter().map(|&s| s.to_string()).collect();

            let res = d4_bgzf::run_tasks(
                files.clone(),
                Some(&sample_names),
                args.clone(),
                progress_bar.clone(),
            )?;

            temp_file_paths.push(io::write_d4_parallel::<PathBuf>(
                &res,
                chroms.clone(),
                None,
            )?);

            let population_name = population_map.get_popname_refs()[idx];
            io::write_bed(
                args.outdir
                    .join(format!("{}_callable_sites.bed", population_name))
                    .as_std_path()
                    .to_path_buf(),
                &res,
            );
        }

        io::merge_d4_files(
            args.outdir
                .join("callable_sites.d4")
                .as_std_path()
                .to_path_buf(),
            temp_file_paths,
            population_map.get_popname_refs(),
        )?;
    } else {
        let res = d4_bgzf::run_tasks(files, None, args.clone(), progress_bar)?;
        io::write_d4_parallel::<PathBuf>(
            &res,
            chroms,
            Some(
                args.outdir
                    .join("callable_sites.d4")
                    .as_std_path()
                    .to_path_buf(),
            ),
        )?;
        io::write_bed(
            args.outdir
                .join("callable_sites.bed")
                .as_std_path()
                .to_path_buf(),
            &res,
        )?;
    }

    Ok(())
}

fn process_merged_d4(
    file: PathBuf,
    args: &LociArgs,
    progress_bar: Option<ProgressBar>,
) -> Result<()> {
    let chroms: Vec<d4::Chrom> = d4_tasks::get_chrom_regions(&file)?
        .iter()
        .map(|r| d4::Chrom {
            name: r.0.to_string(),
            size: r.2.try_into().unwrap(),
        })
        .collect();
    if let Some(population_map) = &args.population_map {
        let mut temp_file_paths = Vec::with_capacity(population_map.num_populations());
        for (idx, samples) in population_map
            .get_samples_per_population()
            .iter()
            .enumerate()
        {
            let sample_names: Vec<String> = samples.iter().map(|&s| s.to_string()).collect();
            let res =
                d4_tasks::run_tasks_on_tracks(file.clone(), Some(sample_names), args.clone())?;

            temp_file_paths.push(io::write_d4_parallel::<PathBuf>(
                &res,
                chroms.clone(),
                None,
            )?);

            let population_name = population_map.get_popname_refs()[idx];
            io::write_bed(
                args.outdir
                    .join(format!("{}_callable_sites.bed", population_name))
                    .as_std_path()
                    .to_path_buf(),
                &res,
            );
        }

        io::merge_d4_files(
            args.outdir
                .join("callable_sites.d4")
                .as_std_path()
                .to_path_buf(),
            temp_file_paths,
            population_map.get_popname_refs(),
        )?;
    } else {
        let res = d4_tasks::run_tasks_on_tracks(file.clone(), None, args.clone())?;
        io::write_d4_parallel::<PathBuf>(
            &res,
            chroms,
            Some(
                args.outdir
                    .join("callable_sites.d4")
                    .as_std_path()
                    .to_path_buf(),
            ),
        )?;
        io::write_bed(
            args.outdir
                .join("callable_sites.bed")
                .as_std_path()
                .to_path_buf(),
            &res,
        )?;
    }

    Ok(())
   
}

// pub enum D4Reader {
//     D4(d4::D4TrackReader),
//     Bgzf(d4_bgzf::BGZID4MatrixReader),
// }

// impl D4Reader {
//     pub fn chrom_regions(
//         &self,
//         thresholds: Thresholds,
//         exclude_chrs: Option<&HashSet<String>>,
//     ) -> Result<Vec<ChromRegion>> {
//         match self {
//             D4Reader::D4(reader) => {
//                 prepare_chrom_regions(reader.chrom_regions(), thresholds, exclude_chrs)
//             }
//             D4Reader::Bgzf(reader) => {
//                 prepare_chrom_regions(reader.chrom_regions(), thresholds, exclude_chrs)
//             }
//         }
//     }

//     pub fn run_tasks(
//         &self,
//         loci_args: &LociArgs,
//         chrom_regions: Vec<ChromRegion>,
//         sample_refs: Option<Vec<String>>,
//         progress_bar: Option<indicatif::ProgressBar>,
//     ) -> Result<Vec<(String, u32, Vec<CallableRegion>)>> {
//         let output_counts = !loci_args.no_counts;
//         let paths = loci_args.read_filelist(None)?;
//         match self {
//             D4Reader::Bgzf(_) => d4_bgzf::run_tasks(
//                 paths,
//                 loci_args.threads,
//                 &chrom_regions,
//                 (loci_args.mean_depth_min, loci_args.mean_depth_max),
//                 loci_args.depth_proportion,
//                 progress_bar,
//             ),
//             D4Reader::D4(_) => {
//                 let tracks = d4_tasks::prepare_tracks_from_file(
//                     loci_args.infile.clone().unwrap(),
//                     sample_refs,
//                 )?;
//                 d4_tasks::run_tasks_on_tracks(
//                     tracks,
//                     chrom_regions,
//                     (loci_args.mean_depth_min, loci_args.mean_depth_max),
//                     loci_args.depth_proportion,
//                     output_counts,
//                 )
//             }
//         }
//     }
// }

// pub enum LociMode {
//     D4,
//     BGZF,
//     GVCF,
// }

// // #[cfg(test)]
// // mod tests {

// //     use std::collections::HashMap;

// //     use super::*;
// //     fn init() {
// //         let _ = env_logger::builder()
// //             .target(env_logger::Target::Stdout)
// //             .filter_level(log::LevelFilter::Trace)
// //             .is_test(true)
// //             .try_init();
// //     }

// //     #[test]
// //     fn test_add_filters_to_chroms_success() {
// //         // Setup input
// //         init();
// //         let d4_chrom_regions = vec![
// //             ("chr1", 0, 1000),
// //             ("chr2", 1000, 2000),
// //             ("chr3", 2000, 3000),
// //         ];
// //         let mut filter_map = HashMap::new();
// //         filter_map.insert("chr1".to_string(), (0.5, 1.5));
// //         filter_map.insert("chr2".to_string(), (0.1, 0.9));
// //         filter_map.insert("chr3".to_string(), (0.2, 1.0));

// //         // Expected output
// //         let expected = vec![
// //             ("chr1".to_string(), 0, 1000, 0.5, 1.5),
// //             ("chr2".to_string(), 1000, 2000, 0.1, 0.9),
// //             ("chr3".to_string(), 2000, 3000, 0.2, 1.0),
// //         ];

// //         // Call function
// //         let result = add_filters_to_chroms(d4_chrom_regions, filter_map).unwrap();

// //         // Assert the result matches expected
// //         assert_eq!(result, expected);
// //     }

// //     #[test]
// //     fn test_add_filters_to_chroms_missing_filter() {
// //         // Setup input with one chromosome missing from the filter map
// //         let d4_chrom_regions = vec![
// //             ("chr1", 0, 1000),
// //             ("chr2", 1000, 2000),
// //             ("chr3", 2000, 3000),
// //         ];
// //         let mut filter_map = HashMap::new();
// //         filter_map.insert("chr1".to_string(), (0.5, 1.5));
// //         filter_map.insert("chr2".to_string(), (0.1, 0.9));
// //         // chr3 is missing in the filter map

// //         // Expected output
// //         let expected = vec![
// //             ("chr1".to_string(), 0, 1000, 0.5, 1.5),
// //             ("chr2".to_string(), 1000, 2000, 0.1, 0.9),
// //             ("chr3".to_string(), 2000, 3000, 0.0, f64::INFINITY), // Default filter for chr3
// //         ];

// //         // Call function
// //         let result = add_filters_to_chroms(d4_chrom_regions, filter_map).unwrap();

// //         // Assert the result matches expected
// //         assert_eq!(result, expected);
// //     }

// //     #[test]
// //     fn test_add_filters_to_chroms_unmatched_filter() {
// //         // Setup input with an extra chromosome in the filter map
// //         let d4_chrom_regions = vec![("chr1", 0, 1000), ("chr2", 1000, 2000)];
// //         let mut filter_map = HashMap::new();
// //         filter_map.insert("chr1".to_string(), (0.5, 1.5));
// //         filter_map.insert("chr2".to_string(), (0.1, 0.9));
// //         filter_map.insert("chrX".to_string(), (0.3, 1.2)); // Unmatched chromosome filter

// //         // Call function, expecting an error
// //         let result = add_filters_to_chroms(d4_chrom_regions, filter_map);

// //         // Assert the function returns an error due to unmatched chromosome in the filter map
// //         assert!(result.is_err());
// //         let err_message = result.unwrap_err().to_string();
// //         assert!(err_message
// //             .contains("Filter provided for chromosome 'chrX' not found in D4 chrom regions."));
// //     }
// // }
