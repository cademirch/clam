mod cli;
mod loci;
mod utils;
use std::{fs::File, path::PathBuf};

use anyhow::{bail, Ok, Result};
use clap::{Arg, Parser, Subcommand};

#[derive(Debug, Parser)]
#[command(name = "clam")]
#[command(about = "Callable Loci and More")]
#[command(version)]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Debug, Subcommand)]
enum Commands {
    Loci(cli::LociArgs),

    Stat {
        file: String,
    },

    #[command(hide = true)]
    Mkdocs,
}
fn main() -> Result<()> {
    env_logger::init();
    let args = Cli::parse();

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

            let bed_out = match loci_args.outfile.extension() {
                Some(ext) if ext == "d4" => false,
                Some(ext) if ext == "bed" => true,
                Some(ext) => bail!("Unrecognized input file extension: {}", ext.to_string()),
                None => bail!("Input file has no extension"),
            };

            if bed_out && populations {
                bail!("To use populations feature, output file must be a d4 file!");
            }

            let d4_reader = if gzipped {
                loci::D4Reader::Bgzf(loci::d4_bgzf::BgzfD4MatrixReader::from_path(
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
                )?;
                let mut temp_file_paths = Vec::with_capacity(population_map.num_populations);

                for samples in population_map.pop_idx_to_sample_names.iter() {
                    let res = d4_reader.run_tasks(
                        &loci_args,
                        chrom_regions.clone(),
                        Some(samples.clone()),
                    )?;
                    temp_file_paths.push(loci::write_d4::<PathBuf>(res, chroms.clone(), None)?);
                }

                loci::merge_d4_files(
                    loci_args.outfile.clone().into_std_path_buf(),
                    temp_file_paths,
                    population_map.get_popname_refs(),
                )?;
            } else {
                let res = d4_reader.run_tasks(&loci_args, chrom_regions.clone(), None)?;
                if bed_out {
                    loci::write_bed(loci_args.outfile.clone(), res, loci_args.output_counts)?;
                } else {
                    loci::write_d4::<PathBuf>(
                        res,
                        chroms,
                        Some(loci_args.outfile.clone().into_std_path_buf()),
                    )?;
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
