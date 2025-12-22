use clam::core::zarr::CallableArrays;
use clap::{Args, Parser, Subcommand};
use color_eyre::eyre::{Context, Ok};
use color_eyre::Result;
use std::collections::HashSet;
use std::path::PathBuf;

/// Population genetics toolkit
#[derive(Parser, Debug)]
#[command(author, version, about = "Population genetics analysis toolkit")]
pub struct Cli {
    #[command(subcommand)]
    pub command: Commands,
}

#[derive(Subcommand, Debug)]
pub enum Commands {
    /// Calculate callable sites from depth statistics
    Loci(LociArgs),
    /// Calculate population genetic statistics from VCF
    Stat(StatArgs),
    /// Collect depth from multiple files into a Zarr store
    Collect(CollectArgs),
}

#[derive(Args, Debug, Clone)]
pub struct SharedOptions {
    /// Number of threads to use for parallel processing
    #[arg(short = 't', long = "threads", default_value_t = 1)]
    pub threads: usize,

    /// Path to file that defines populations. Tab separated: sample\tpopulation_name
    #[arg(short = 'p', long = "population-file")]
    pub population_file: Option<PathBuf>,

    /// Comma separated list of chromosomes to exclude
    #[arg(short = 'x', long = "exclude", value_delimiter = ',', num_args = 1.., conflicts_with = "exclude_file")]
    pub exclude: Option<Vec<String>>,

    /// Path to file with chromosomes to exclude, one per line
    #[arg(long = "exclude-file")]
    pub exclude_file: Option<PathBuf>,

    /// Comma separated list of chromosomes to include (restrict analysis to)
    #[arg(short = 'i', long = "include", value_delimiter = ',', num_args = 1.., conflicts_with = "include_file")]
    pub include: Option<Vec<String>>,

    /// Path to file with chromosomes to include, one per line
    #[arg(long = "include-file")]
    pub include_file: Option<PathBuf>,
}

impl SharedOptions {
    /// Initialize thread pool
    pub fn initialize_threading(&self) -> Result<()> {
        rayon::ThreadPoolBuilder::new()
            .num_threads(self.threads)
            .build_global()?;
        Ok(())
    }

    /// Get excluded chromosomes as a HashSet
    pub fn get_excluded_chromosomes(&self) -> Result<Option<HashSet<String>>> {
        if let Some(ref exclude_list) = self.exclude {
            return Ok(Some(exclude_list.iter().cloned().collect()));
        }

        if let Some(ref exclude_file) = self.exclude_file {
            let content = std::fs::read_to_string(exclude_file)?;
            let chroms = content
                .lines()
                .map(|s| s.trim().to_string())
                .filter(|s| !s.is_empty())
                .collect();
            return Ok(Some(chroms));
        }

        Ok(None)
    }

    /// Get included chromosomes as a HashSet
    pub fn get_included_chromosomes(&self) -> Result<Option<HashSet<String>>> {
        if let Some(ref include_list) = self.include {
            return Ok(Some(include_list.iter().cloned().collect()));
        }

        if let Some(ref include_file) = self.include_file {
            let content = std::fs::read_to_string(include_file)?;
            let chroms = content
                .lines()
                .map(|s| s.trim().to_string())
                .filter(|s| !s.is_empty())
                .collect();
            return Ok(Some(chroms));
        }

        Ok(None)
    }
}

/// Calculate callable sites from depth statistics
#[derive(Args, Debug)]
pub struct LociArgs {
    /// Input depth files
    #[arg(required = true)]
    pub input: Vec<PathBuf>,

    /// Output path for callable sites zarr array
    #[arg(short = 'o', long = "output", required = true)]
    pub output: PathBuf,

    /// Minimum depth to consider a site callable for each individual
    #[arg(short = 'm', long = "min-depth", default_value_t = 0.0)]
    pub min_depth: f64,

    /// Maximum depth to consider a site callable for each individual
    #[arg(short = 'M', long = "max-depth", default_value_t = f64::INFINITY)]
    pub max_depth: f64,

    /// Proportion of samples that must pass thresholds at a site
    #[arg(short = 'd', long = "min-proportion", default_value_t = 0.0)]
    pub min_proportion: f64,

    /// Minimum mean depth across all samples required at a site
    #[arg(long = "min-mean-depth", default_value_t = 0.0)]
    pub mean_depth_min: f64,

    /// Maximum mean depth across all samples allowed at a site
    #[arg(long = "max-mean-depth", default_value_t = f64::INFINITY)]
    pub mean_depth_max: f64,

    /// Chunk size for processing (base pairs)
    #[arg(long = "chunk-size", default_value_t = 1_000_000)]
    pub chunk_size: u64,

    /// Output per-sample masks instead of per-population counts
    #[arg(long = "per-sample")]
    pub per_sample: bool,

    /// Minimum gq to count depth (GVCF input only)
    #[arg(long = "min-gq")]
    pub min_gq: Option<isize>,

    /// Shared options
    #[command(flatten)]
    pub shared: SharedOptions,
}

#[derive(Args, Debug, Clone)]
pub struct CollectArgs {
    #[arg(required = true)]
    pub input: Vec<PathBuf>,

    /// Output path for callable sites zarr array
    #[arg(short = 'o', long = "output", required = true)]
    pub output: PathBuf,

    /// Chunk size for processing (base pairs)
    #[arg(long = "chunk-size", default_value_t = 1_000_000)]
    pub chunk_size: u64,

    /// Minimum gq to count depth (GVCF input only)
    #[arg(long = "min-gq")]
    pub min_gq: Option<isize>,

    /// Shared options
    #[command(flatten)]
    pub shared: SharedOptions,
}

/// Calculate population genetic statistics from VCF
#[derive(Args, Debug)]
pub struct StatArgs {
    /// Path to input VCF file (bgzipped and indexed)
    #[arg(required = true)]
    pub vcf: PathBuf,

    /// Output directory for statistics files
    #[arg(short = 'o', long = "outdir", required = true)]
    pub outdir: PathBuf,

    /// Path to callable sites zarr array (from clam loci)
    #[arg(short = 'c', long = "callable")]
    pub callable: Option<PathBuf>,

    /// Path to ROH regions BED file (sample name in 4th column)
    #[arg(short = 'r', long = "roh")]
    pub roh: Option<PathBuf>,

    /// Window size in base pairs
    #[arg(short = 'w', long = "window-size", conflicts_with = "regions_file")]
    pub window_size: Option<usize>,

    /// BED file specifying regions to calculate statistics for
    #[arg(long = "regions-file", conflicts_with = "window_size")]
    pub regions_file: Option<PathBuf>,

    /// Chunk size for parallel processing (base pairs)
    #[arg(long = "chunk-size", default_value_t = 1_000_000)]
    pub chunk_size: u64,

    /// Shared options
    #[command(flatten)]
    pub shared: SharedOptions,
}

impl CollectArgs {
    pub fn run(self) -> Result<()> {
        use clam::collect::run_collect;
        run_collect(self.input, self.output, self.chunk_size, self.min_gq)?;
        Ok(())
    }
}

impl LociArgs {
    pub fn run(self) -> Result<()> {
        use clam::core::population::PopulationMap;
        use clam::core::zarr::is_zarr_path;
        use clam::loci::{run_loci, run_loci_zarr, ThresholdConfig};

        self.shared.initialize_threading()?;

        let pop_map = if let Some(ref pop_file) = self.shared.population_file {
            Some(PopulationMap::from_file(pop_file)?)
        } else {
            None
        };

        let thresholds = ThresholdConfig {
            min_depth: self.min_depth,
            max_depth: self.max_depth,
            min_proportion: self.min_proportion,
            mean_depth_range: (self.mean_depth_min, self.mean_depth_max),
        };

        if self.input.len() == 1 && is_zarr_path(&self.input[0]) {
            let zarr = &self.input[0];
            run_loci_zarr(
                zarr.to_path_buf(),
                self.output,
                pop_map,
                thresholds,
                self.per_sample,
            )
        } else {
            run_loci(
                self.input,
                self.output,
                pop_map,
                thresholds,
                self.chunk_size,
                self.per_sample,
                self.min_gq,
            )
        }
    }
}

impl StatArgs {
    pub fn run(self) -> Result<()> {
        use clam::core::zarr::is_zarr_path;
        use clam::stat::config::StatConfig;
        use clam::stat::run_stat;
        use clam::stat::utils::read_bed_regions;
        use clam::stat::windows::WindowStrategy;

        self.shared.initialize_threading()?;

        std::fs::create_dir_all(&self.outdir).wrap_err(format!(
            "Failed to create output directory: {}",
            self.outdir.display()
        ))?;

        let window_strategy = if let Some(size) = self.window_size {
            WindowStrategy::FixedSize(size)
        } else if let Some(ref regions_file) = self.regions_file {
            let regions = read_bed_regions(regions_file)?;
            WindowStrategy::Regions(regions)
        } else {
            return Err(color_eyre::eyre::eyre!(
                "Must specify either --window-size or --regions-file"
            ));
        };
        let exclude = self.shared.get_excluded_chromosomes()?;
        let include = self.shared.get_included_chromosomes()?;

        let mut chunk_size = self.chunk_size;
        if let Some(callable_path) = self.callable.as_ref() {
            if is_zarr_path(&callable_path) {
                let callable_zarr = CallableArrays::open(&callable_path)?;
                chunk_size = callable_zarr.chunk_size();
            }
        }

        let config = StatConfig::new(
            self.vcf,
            self.shared.population_file,
            self.callable,
            self.roh,
            window_strategy,
            chunk_size,
            exclude,
            include,
        )?;

        run_stat(config, &self.outdir)
    }
}

// Main entry point
pub fn main() -> Result<()> {
    color_eyre::install()?;
    use env_logger::Env;

    env_logger::Builder::from_env(Env::default().default_filter_or("info")).init();
    let cli = Cli::parse();

    match cli.command {
        Commands::Loci(args) => args.run(),
        Commands::Stat(args) => args.run(),
        Commands::Collect(args) => args.run(),
    }
}
