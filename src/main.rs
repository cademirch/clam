use clam::{loci, stat};
use anyhow::Result;
use clap::{Parser, Subcommand};
use indicatif::{MultiProgress, ProgressBar, ProgressStyle};
use indicatif_log_bridge::LogWrapper;
use log::info;

#[derive(Debug, Parser)]
#[command(name = "clam")]
#[command(about = "Callable Loci and More")]
#[command(version)]
struct Cli {
    #[command(subcommand)]
    command: Commands,

    /// Increase verbosity (-v, -vv for more verbosity)
    #[arg(short, long, action = clap::ArgAction::Count)]
    verbose: u8,

    /// Suppress output (overrides verbosity)
    #[arg(short, long, action)]
    quiet: bool,
}
#[derive(Debug, Subcommand)]
enum Commands {
    Loci(loci::LociArgs),

    Stat(stat::StatArgs),

    #[command(hide = true)]
    Mkdocs,
}
fn main() -> Result<()> {
    let args = Cli::parse();

    // Set up logging based on verbosity flags
    let log_level = if args.quiet {
        log::LevelFilter::Off
    } else {
        match args.verbose {
            0 => log::LevelFilter::Info,  // Default
            1 => log::LevelFilter::Debug, // -v
            2 => log::LevelFilter::Trace, // -vv
            _ => log::LevelFilter::Trace, // Any more `v`s
        }
    };

    let logger = env_logger::builder().filter_level(log_level).build();
    let multi = MultiProgress::new();
    LogWrapper::new(multi.clone(), logger).try_init().unwrap();
    log::set_max_level(log_level);
    let progress_bar = if args.quiet {
        None
    } else {
        let bar = ProgressBar::no_length().with_style(
            ProgressStyle::with_template(
                "{spinner:.green} [{wide_bar:.cyan/blue}] {percent}% done.",
            )
            .unwrap()
            .progress_chars("#>-"),
        );
        Some(multi.add(bar))
    };

    // Capture command-line arguments for logging
    let str_args: Vec<String> = std::env::args().collect();
    let command_line = str_args.join(" ");
    let version = env!("CARGO_PKG_VERSION");
    info!(
        "clam version: {} arguments supplied: {}",
        version, command_line
    );
    match args.command {
        Commands::Mkdocs => {
            clap_markdown::print_help_markdown::<Cli>();
            return Ok(());
        }
        Commands::Loci(loci_args) => {
            loci::process(loci_args, progress_bar)?;

            Ok(())
        }
        Commands::Stat(stat_args) => {
            stat::run_stat(stat_args, progress_bar)?;
            Ok(())
        }
    }
}
