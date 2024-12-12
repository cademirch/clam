mod loci;
mod stat;
mod utils;

use std::fs::create_dir_all;
use std::path::PathBuf;

use anyhow::{Context, Result};
use clap::{Parser, Subcommand};
use indicatif::{MultiProgress, ProgressBar, ProgressStyle};
use indicatif_log_bridge::LogWrapper;
use log::{info, warn};

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
            let _thread_pool = rayon::ThreadPoolBuilder::new()
                .num_threads(loci_args.threads)
                .build_global()?;

            let thresholds = loci_args.threshold_file.as_ref().map_or_else(
                || loci::Thresholds::Fixed((loci_args.min_depth, loci_args.max_depth)),
                |file| {
                    loci::read_threshold_file(file)
                        .map(loci::Thresholds::PerChromosome)
                        .unwrap_or_else(|_| {
                            loci::Thresholds::Fixed((loci_args.min_depth, loci_args.max_depth))
                        })
                },
            );

            let exclude_chrs =
                utils::get_exclude_chromosomes(&loci_args.exclude, &loci_args.exclude_file)?;
            let gzipped = loci_args.infile.extension().unwrap() == "gz";
            let populations = loci_args.population_file.is_some();

            let bed_out = loci_args.no_counts;
            let outfile = loci_args.resolve_output_file();

            if let Some(parent) = outfile.parent() {
                if !parent.exists() {
                    create_dir_all(parent).with_context(|| {
                        format!("Failed to create parent directories for {:?}", outfile)
                    })?;
                }
            }
            if !loci_args.no_counts
                && (loci_args.mean_depth_min > 0.0 || loci_args.depth_proportion > 0.0)
            {
                warn!("Mean depth proportion settings ignored because we are reporting counts.");
            }

            let d4_reader = if gzipped {
                loci::D4Reader::Bgzf(loci::d4_bgzf::BGZID4MatrixReader::from_merged(
                    loci_args.infile.clone(),
                    None,
                )?)
            } else {
                loci::D4Reader::D4(d4::D4TrackReader::open_first_track(
                    loci_args.infile.clone().into_std_path_buf(),
                )?)
            };

            let chrom_regions = d4_reader.chrom_regions(thresholds.clone(), Some(&exclude_chrs))?;
            let chroms: Vec<d4::Chrom> = chrom_regions
                .iter()
                .map(|r| d4::Chrom {
                    name: r.chr.to_string(),
                    size: r.end.try_into().unwrap(),
                })
                .collect();

            if populations {
                let population_map = utils::PopulationMapping::from_path(
                    loci_args.population_file.as_ref().unwrap(),
                    None,
                )?;
                let mut temp_file_paths = Vec::with_capacity(population_map.num_populations);

                for samples in population_map.pop_idx_to_sample_names.iter() {
                    let res = d4_reader.run_tasks(
                        &loci_args,
                        chrom_regions.clone(),
                        Some(samples.clone()),
                    )?;
                    temp_file_paths.push(loci::write_d4_parallel::<PathBuf>(
                        res,
                        chroms.clone(),
                        None,
                    )?);
                }

                loci::merge_d4_files(
                    outfile.clone().into_std_path_buf(),
                    temp_file_paths,
                    population_map.get_popname_refs().unwrap(),
                )?;
            } else {
                let res = d4_reader.run_tasks(&loci_args, chrom_regions.clone(), None)?;
                if bed_out {
                    loci::write_bed(outfile.clone(), res)?;
                } else {
                    loci::write_d4_parallel::<PathBuf>(
                        res,
                        chroms,
                        Some(outfile.clone().into_std_path_buf()),
                    )?;
                }
            }

            Ok(())
        }
        Commands::Stat(stat_args) => {
            stat::run_stat(stat_args, progress_bar)?;
            Ok(())
        }
    }
}
