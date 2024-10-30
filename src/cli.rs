use clap::{Parser, Subcommand, ArgGroup};
use camino::Utf8PathBuf;

#[derive(Debug, Parser)]
#[command(name = "clam")]
#[command(about = "Callable Loci and More")]
pub struct Cli {
    #[command(subcommand)]
    pub command: Commands,
}

#[derive(Debug, Subcommand)]
pub enum Commands {
    /// Loci calculation subcommand with detailed depth-based arguments
    Loci(LociArgs),

    /// Stat calculation subcommand
    Stat {
        /// Path to input file
        file: String,
    },
}

#[derive(Parser, Debug, Clone)]
#[command(author, version, about = "Calculates callable sites from depth statistics.")]
#[command(
    group(
        ArgGroup::new("input")
            .required(true)
            .args(&["infile"])
    ),
    group(
        ArgGroup::new("output")
            .required(true)
            .args(&["outprefix"])
    ),
    group(
        ArgGroup::new("global_thresholds")
            .args(&["min_depth", "max_depth"])
    ),
    group(
        ArgGroup::new("threshold_file")
            .args(&["threshold_file"])
            .conflicts_with("global_thresholds")
    ),
    group(
        ArgGroup::new("populations")
            .args(&["population_file"])
            .conflicts_with("nopops")
    ),
    group(
        ArgGroup::new("nopops")
            .args(&["depth_proportion", "mean_depth_min", "mean_depth_max"])
    ),
    group(
        ArgGroup::new("exclude_group")
        .args(&["exclude", "exclude-file"])
        .required(false)
        .multiple(false)
    )
)]
pub struct LociArgs {
    /// Path to input D4 file
    pub infile: Utf8PathBuf,
    /// Prefix for output
    pub outprefix: Utf8PathBuf,
    /// Minimum depth to consider site callable per individual
    #[arg(short = 'm', long = "min-depth", default_value_t = 0.0)]
    pub min_depth: f64,
    /// Maximum depth to consider site callable per individual
    #[arg(short = 'M', long = "max-depth", default_value_t = f64::INFINITY)]
    pub max_depth: f64,
    /// Proportion of samples passing thresholds at site to consider callable
    #[arg(short = 'd', long = "depth-proportion", default_value_t = 1.0)]
    pub depth_proportion: f64,
    /// Minimum mean depth across all samples at site to consider callable
    #[arg(short = 'u', long = "min-mean-depth", default_value_t = 0.0)]
    pub mean_depth_min: f64,
    /// Maximum mean depth across all samples at site to consider callable
    #[arg(short = 'U', long = "max-mean-depth", default_value_t = f64::INFINITY)]
    pub mean_depth_max: f64,
    /// Output number of individuals callable at site. EXPERIMENTAL v0.1.0
    #[arg(short = 'c', long = "output-counts", default_value_t = false)]
    pub output_counts: bool,
    /// Number of threads to use
    #[arg(short = 't', long = "threads", default_value_t = 1)]
    pub threads: usize,
    /// Path to file that defines populations. Tab separated: sample, population_name
    #[arg(short = 'p', long = "population-file")]
    pub population_file: Option<Utf8PathBuf>,
    /// Path to file that defines per-chromosome individual level thresholds. Tab separated: chrom, min, max
    #[arg(short = 't', long = "thresholds-file")]
    pub threshold_file: Option<Utf8PathBuf>,
    /// Comma separated list of chromosomes to exclude
    #[arg(short = 'x', value_delimiter = ',', num_args = 1..)]
    pub exclude: Option<Vec<String>>,
    /// Path to file with chromosomes to exclude, one per line
    #[arg(long = "exclude-file")]
    pub exclude_file: Option<Utf8PathBuf>,
}
