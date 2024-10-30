mod utils;
mod d4_tasks;
mod d4_bgzf;

use clap::{Args, Parser, Subcommand, ArgGroup};
use d4_bgzf::BgzfD4MatrixReader;
use camino::Utf8PathBuf;
use rayon::ThreadPoolBuilder;

#[derive(Debug, Parser)]
#[command(name = "clam")]
#[command(about = "Callable Loci and More")]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}



#[derive(Debug, Subcommand)]
enum Commands {
    /// Loci calculation subcommand with detailed depth-based arguments
    Loci(LociArgs),

    /// Stat calculation subcommand
    Stat {
        /// Path to input file
        file: String,
    },
}

#[derive(Parser, Debug)]
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
    )
)]
struct LociArgs {
    /// Path to input D4 file
    infile: Utf8PathBuf,

    /// Prefix for output
    outprefix: Utf8PathBuf,

    /// Minimum depth to consider site callable per individual
    #[arg(short = 'm', long = "min-depth", default_value_t = 0.0)]
    min_depth: f64,

    /// Maximum depth to consider site callable per individual
    #[arg(short = 'M', long = "max-depth", default_value_t = f64::INFINITY)]
    max_depth: f64,

    /// Proportion of samples passing thresholds at site to consider callable
    #[arg(short = 'd', long = "depth-proportion", default_value_t = 1.0)]
    depth_proportion: f64,

    /// Minimum mean depth across all samples at site to consider callable
    #[arg(short = 'u', long = "min-mean-depth", default_value_t = 0.0)]
    mean_depth_min: f64,

    /// Maximum mean depth across all samples at site to consider callable
    #[arg(short = 'U', long = "max-mean-depth", default_value_t = f64::INFINITY)]
    mean_depth_max: f64,

    /// Output number of individuals callable at site. EXPERIMENTAL v0.1.0
    #[arg(short = 'c', long = "output-counts", default_value_t = false)]
    output_counts: bool,

    /// Number of threads to use
    #[arg(short = 't', long = "threads", default_value_t = 1)]
    threads: usize,

    /// Path to file that defines populations. Tab separated: sample, population_name
    #[arg(short = 'p', long = "population-file")]
    population_file: Option<Utf8PathBuf>,

    /// Path to file that defines per-chromosome individual level thresholds. Tab separated: chrom, min, max
    #[arg(short = 't', long = "thresholds-file")]
    threshold_file: Option<Utf8PathBuf>,

    /// Comma separated list of chromosomes to exclude
    #[arg(short = 'x', value_delimiter = ',', num_args = 1..)]
    exclude: Option<Vec<String>>,

    /// Path to file with chromosomes to exclude, one per line
    #[arg(long = "exclude-file")]
    exclude_file: Option<Utf8PathBuf>,
}

fn main() {
    let args = Cli::parse();

    match args.command {
        Commands::Loci(loci_args) => {
            // Initialize the thread pool
            let _tpool = ThreadPoolBuilder::new()
                .num_threads(loci_args.threads)
                .build_global()
                .expect("Failed to build thread pool");

            // Create the matrix reader
            let mut matrix_rdr = BgzfD4MatrixReader::from_path(&loci_args.infile, None)
                .expect("Failed to create BgzfD4MatrixReader");

            // Generate a completely independent collection of regions
            let regions: Vec<(String, u32, u32)> = matrix_rdr
                .chrom_regions()
                .into_iter()
                .map(|(chrom, start, end)| (chrom.to_string(), start, end))
                .collect();

            // Now, we can safely iterate with mutable access to `matrix_rdr`
            for (i, (chrom, begin, end)) in regions.iter().enumerate() {
                let callable = matrix_rdr
                    .get_callable_regions(
                        chrom,
                        *begin,
                        *end,
                        (loci_args.min_depth, loci_args.max_depth),
                        (loci_args.mean_depth_min, loci_args.mean_depth_max),
                        loci_args.depth_proportion,
                        true,
                    )
                    .expect("Failed to get callable regions");

                // Process `callable` as needed
                println!("Region {} done", i);
            }
        }
        Commands::Stat { file } => {
            // Placeholder for the Stat command functionality
            todo!();
        }
    }
}