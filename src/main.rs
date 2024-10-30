mod cli;
mod loci;
mod utils;
use std::fs::File;

use anyhow::Result;
use clap::{Parser, Subcommand};

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
    Loci(cli::LociArgs),

    /// Stat calculation subcommand
    Stat {
        /// Path to input file
        file: String,
    },
}

fn main() -> Result<()> {
    env_logger::init();
    let args = Cli::parse();

    match args.command {
        Commands::Loci(loci_args) => {
            let thresholds = if let Some(threshold_file) = &loci_args.threshold_file {
                let chrom_thresholds = loci::read_threshold_file(threshold_file)?;
                loci::Thresholds::PerChromosome(chrom_thresholds)
            } else {
                loci::Thresholds::Fixed((loci_args.min_depth, loci_args.max_depth))
            };
            let exclude_chrs =
                utils::get_exclude_chromosomes(&loci_args.exclude, &loci_args.exclude_file);
            let gzipped = loci_args.infile.extension().unwrap() == "gz";
            let populations = loci_args.population_file.is_some();

            if gzipped {
                let mut rdr = loci::d4_bgzf::BgzfD4MatrixReader::from_path(
                    loci_args.infile.to_path_buf(),
                    None,
                )?;
                let chrom_regions = loci::prepare_chrom_regions(
                    rdr.chrom_regions(),
                    thresholds,
                    Some(&exclude_chrs?),
                )?;
                drop(rdr);
                if populations {
                    let population_map =
                        utils::PopulationMapping::from_path(loci_args.population_file.unwrap())?;
                    for (pop_idx, samples) in
                        population_map.pop_idx_to_sample_names.iter().enumerate()
                    {
                        loci::d4_bgzf::run_bgzf_tasks(
                            loci_args.infile.to_path_buf(),
                            loci_args.threads,
                            samples.to_vec(),
                            chrom_regions.clone(),
                            (loci_args.mean_depth_min, loci_args.mean_depth_max),
                            loci_args.depth_proportion,
                            loci_args.output_counts,
                        )?;
                    }
                }
            }

            Ok(())
        }
        Commands::Stat { file } => {
            // Placeholder for the Stat command functionality
            todo!();
        }
    }
}
