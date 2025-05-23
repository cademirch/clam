pub mod counts;
pub mod d4_bgzf;
pub mod d4_tasks;
pub mod gvcf;
pub mod io;
pub mod regions;
pub mod thresholds;

use std::collections::{HashSet, VecDeque};
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::PathBuf;

use anyhow::{Context, Result};
use camino::Utf8PathBuf;
use clap::Parser;
use indicatif::ProgressBar;


#[derive(Debug, Clone, PartialEq)]
pub enum InputMode {
    MergedD4(PathBuf),            // Single merged D4 file
    MultiD4(VecDeque<PathBuf>),   // Multiple per-sample D4 files
    MultiGVCF(VecDeque<PathBuf>), // Multiple per-sample GVCF files
}

#[derive(Parser, Debug, Clone)]
#[command(
    author,
    version,
    about = "Calculate callable sites from depth statistics.",
    // Add some examples of common use cases
    after_help = "EXAMPLES:
    # Basic usage with positional input files
    clam loci -o output_dir input1.d4 input2.d4
    
    # Using a file list instead of positional arguments
    clam loci -f filelist.txt -o output_dir
    
    # Set custom depth thresholds
    clam loci -o output_dir -m 10 -M 100 input1.d4 input2.d4"
)]
pub struct LociArgs {
    // === Input Options ===
    /// Input files (D4 format by default). Specify one or more files directly.
    #[arg(required_unless_present = "filelist")]
    pub input: Vec<Utf8PathBuf>,

    /// Path to file containing list of input files, one per line.
    /// Use this instead of positional input arguments for many files.
    #[arg(short = 'f', long = "filelist")]
    pub filelist: Option<Utf8PathBuf>,

    /// Use GVCF format instead of default D4 format for input files.
    #[arg(long = "gvcf")]
    pub gvcf: bool,

    /// Input is a merged D4 file (single file containing multiple samples).
    #[arg(long = "merged", conflicts_with = "gvcf")]
    pub merged: bool,

    // === Output Options ===
    /// Output directory for results (required).
    #[arg(short = 'o', required(true))]
    pub outdir: Utf8PathBuf,

    /// Write additional BED file. Note: This can be slow for large datasets.
    #[arg(long = "bed")]
    pub write_bed: bool,

    // === Sample-level Threshold Options ===
    /// Minimum depth to consider a site callable for each individual.
    #[arg(
        short = 'm',
        long = "min-depth",
        default_value_t = 0.0,
        conflicts_with("threshold_file"),
        help_heading = "Sample-level Thresholds"
    )]
    pub min_depth: f64,

    /// Maximum depth to consider a site callable for each individual.
    #[arg(
        short = 'M', 
        long = "max-depth", 
        default_value_t = f64::INFINITY, 
        conflicts_with("threshold_file"),
        help_heading = "Sample-level Thresholds"
    )]
    pub max_depth: f64,

    /// Custom thresholds per chromosome. Tab-separated file: chrom, min, max
    #[arg(
        long = "thresholds-file",
        help_heading = "Sample-level Thresholds"
    )]
    pub threshold_file: Option<Utf8PathBuf>,

    // === Population-level Threshold Options ===
    /// Proportion of samples that must pass thresholds at a site to consider it callable.
    /// Value between 0.0 and 1.0.
    #[arg(
        short = 'd',
        long = "depth-proportion",
        default_value_t = 0.0,
        conflicts_with("population_file"),
        help_heading = "Population-level Thresholds"
    )]
    pub depth_proportion: f64,

    /// Minimum mean depth across all samples required at a site to consider it callable.
    #[arg(
        short = 'u',
        long = "min-mean-depth",
        default_value_t = 0.0,
        conflicts_with("population_file"),
        help_heading = "Population-level Thresholds"
    )]
    pub mean_depth_min: f64,

    /// Maximum mean depth across all samples allowed at a site to consider it callable.
    #[arg(
        short = 'U', 
        long = "max-mean-depth", 
        default_value_t = f64::INFINITY, 
        conflicts_with("population_file"),
        help_heading = "Population-level Thresholds"
    )]
    pub mean_depth_max: f64,

    /// Path to file that defines populations. Tab separated: sample, population_name
    #[arg(
        short = 'p', 
        long = "population-file",
        help_heading = "Population-level Thresholds"
    )]
    pub population_file: Option<Utf8PathBuf>,

    // === Chromosome Filtering Options ===
    /// Comma separated list of chromosomes to exclude.
    /// Example: --exclude chr1,chr2,chrX
    #[arg(
        short = 'x', 
        value_delimiter = ',', 
        num_args = 1.., 
        conflicts_with("exclude_file"),
        help_heading = "Chromosome Filtering"
    )]
    pub exclude: Option<Vec<String>>,

    /// Path to file with chromosomes to exclude, one per line.
    #[arg(
        long = "exclude-file", 
        conflicts_with("exclude"),
        help_heading = "Chromosome Filtering"
    )]
    pub exclude_file: Option<Utf8PathBuf>,

    /// Comma separated list of chromosomes to include (restrict analysis to).
    /// Example: --include chr1,chr2,chr3
    #[arg(
        short = 'i', 
        value_delimiter = ',', 
        num_args = 1.., 
        conflicts_with("include_file"),
        help_heading = "Chromosome Filtering"
    )]
    pub include: Option<Vec<String>>,

    /// Path to file with chromosomes to include, one per line.
    #[arg(
        long = "include-file", 
        conflicts_with("include"),
        help_heading = "Chromosome Filtering"
    )]
    pub include_file: Option<Utf8PathBuf>,

    // === Performance Options ===
    /// Number of threads to use for parallel processing.
    #[arg(
        short = 't', 
        long = "threads", 
        default_value_t = 1,
        help_heading = "Performance"
    )]
    pub threads: usize,

    // Fields that get initialized (no change needed)
    #[arg(skip)]
    pub exclude_chrs: Option<HashSet<String>>,

    #[arg(skip)]
    pub include_chrs: Option<HashSet<String>>,

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

        self.include_chrs =
            super::utils::get_exclude_chromosomes(&self.include, &self.include_file)?;

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

            let res = gvcf::run_tasks(
                files.clone(),
                Some(&sample_names),
                &args,
                progress_bar.clone(),
            )?;

            temp_file_paths.push(io::write_d4_parallel::<PathBuf>(
                &res,
                chroms.clone(),
                None,
            )?);

            let population_name = population_map.get_popname_refs()[idx];
            if args.write_bed {
                let _ = io::write_bed(
                    args.outdir
                        .join(format!("{}_callable_sites.bed", population_name))
                        .as_std_path()
                        .to_path_buf(),
                    &res,
                );
            };
        }

        let _ = io::merge_d4_files(
            args.outdir
                .join("callable_sites.d4")
                .as_std_path()
                .to_path_buf(),
            temp_file_paths,
            population_map.get_popname_refs(),
        )?;
    } else {
        let res = gvcf::run_tasks(files, None, &args, progress_bar)?;
        let _ = io::write_d4_parallel::<PathBuf>(
            &res,
            chroms,
            Some(
                args.outdir
                    .join("callable_sites.d4")
                    .as_std_path()
                    .to_path_buf(),
            ),
        )?;
        if args.write_bed {
            let _ = io::write_bed(
                args.outdir
                    .join("callable_sites.bed")
                    .as_std_path()
                    .to_path_buf(),
                &res,
            )?;
        };
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
            if args.write_bed {
                let _ = io::write_bed(
                    args.outdir
                        .join(format!("{}_callable_sites.bed", population_name))
                        .as_std_path()
                        .to_path_buf(),
                    &res,
                );
            };
        }

        let _ = io::merge_d4_files(
            args.outdir
                .join("callable_sites.d4")
                .as_std_path()
                .to_path_buf(),
            temp_file_paths,
            population_map.get_popname_refs(),
        )?;
    } else {
        let res = d4_bgzf::run_tasks(files, None, args.clone(), progress_bar)?;
        let _ = io::write_d4_parallel::<PathBuf>(
            &res,
            chroms,
            Some(
                args.outdir
                    .join("callable_sites.d4")
                    .as_std_path()
                    .to_path_buf(),
            ),
        )?;
        if args.write_bed {
            let _ = io::write_bed(
                args.outdir
                    .join("callable_sites.bed")
                    .as_std_path()
                    .to_path_buf(),
                &res,
            )?;
        };
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
            if args.write_bed {
                let _ = io::write_bed(
                    args.outdir
                        .join(format!("{}_callable_sites.bed", population_name))
                        .as_std_path()
                        .to_path_buf(),
                    &res,
                );
            }
        }

        let _ = io::merge_d4_files(
            args.outdir
                .join("callable_sites.d4")
                .as_std_path()
                .to_path_buf(),
            temp_file_paths,
            population_map.get_popname_refs(),
        )?;
    } else {
        let res = d4_tasks::run_tasks_on_tracks(file.clone(), None, args.clone())?;
        let _ = io::write_d4_parallel::<PathBuf>(
            &res,
            chroms,
            Some(
                args.outdir
                    .join("callable_sites.d4")
                    .as_std_path()
                    .to_path_buf(),
            ),
        )?;
        let _ = io::write_bed(
            args.outdir
                .join("callable_sites.bed")
                .as_std_path()
                .to_path_buf(),
            &res,
        )?;
    }

    Ok(())
}
