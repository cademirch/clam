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
