mod cli;
mod loci;
mod utils;
use std::{fs::File, path::PathBuf};

use anyhow::{bail, Ok, Result};
use clap::{Arg, Parser, Subcommand};

#[derive(Debug, Parser)]
#[command(name = "clam",)]
#[command(about = "Callable Loci and More")]
#[command(version)]
struct Cli {
    
    #[command(subcommand)]
    command: Commands,
    
}

#[derive(Debug, Subcommand)]
enum Commands {
    Loci(cli::LociArgs),

    Stat { file: String },
    
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
            let bed_out = match loci_args.outfile.extension() {
                Some(ext) if ext == "d4" => false,
                Some(ext) if ext == "bed" => true,
                Some(ext) => bail!("Unrecognized input file extension: {}", ext.to_string()),
                None => bail!("Input file has no extension"),
            };

            if bed_out && populations {
                bail!("To use populations feature, output file must be a d4 file!")
            }

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
                let chroms: Vec<d4::Chrom> = chrom_regions
                    .clone()
                    .into_iter()
                    .map(|r| d4::Chrom {
                        name: r.chr.to_string(),
                        size: r.end.try_into().unwrap(),
                    })
                    .collect();
                drop(rdr);
                if populations {
                    let population_map =
                        utils::PopulationMapping::from_path(loci_args.population_file.unwrap())?;
                    let mut temp_file_paths = Vec::with_capacity(population_map.num_populations);
                    for (_, samples) in population_map.pop_idx_to_sample_names.iter().enumerate() {
                        let res = loci::d4_bgzf::run_bgzf_tasks(
                            loci_args.infile.to_path_buf(),
                            loci_args.threads,
                            Some(samples.to_vec()),
                            chrom_regions.clone(),
                            (loci_args.mean_depth_min, loci_args.mean_depth_max),
                            loci_args.depth_proportion,
                            loci_args.output_counts,
                        )?;
                        let fp = loci::write_d4::<PathBuf>(res, chroms.clone(), None)?;
                        temp_file_paths.push(fp);
                    }

                    loci::merge_d4_files(
                        loci_args.outfile.to_path_buf().into_std_path_buf(),
                        temp_file_paths,
                        population_map.get_popname_refs(),
                    )?;
                } else {
                    let res = loci::d4_bgzf::run_bgzf_tasks(
                        loci_args.infile.to_path_buf(),
                        loci_args.threads,
                        None,
                        chrom_regions.clone(),
                        (loci_args.mean_depth_min, loci_args.mean_depth_max),
                        loci_args.depth_proportion,
                        loci_args.output_counts,
                    )?;

                    if bed_out {
                        loci::write_bed(loci_args.outfile, res, loci_args.output_counts)?;
                    } else {
                        loci::write_d4::<PathBuf>(
                            res,
                            chroms.clone(),
                            Some(loci_args.outfile.into_std_path_buf()),
                        )?;
                    }
                }
            } else {
                let rdr: d4::D4TrackReader = d4::D4TrackReader::open_first_track(
                    loci_args.infile.clone().into_std_path_buf(),
                )?;
                let chrom_regions = loci::prepare_chrom_regions(
                    rdr.chrom_regions(),
                    thresholds,
                    Some(&exclude_chrs?),
                )?;
                let chroms: Vec<d4::Chrom> = chrom_regions
                    .clone()
                    .into_iter()
                    .map(|r| d4::Chrom {
                        name: r.chr.to_string(),
                        size: r.end.try_into().unwrap(),
                    })
                    .collect();
                drop(rdr);
                if populations {
                    let population_map =
                        utils::PopulationMapping::from_path(loci_args.population_file.unwrap())?;
                    let mut temp_file_paths = Vec::with_capacity(population_map.num_populations);

                    for (_, samples) in population_map.pop_idx_to_sample_names.iter().enumerate() {
                        let sample_refs: Vec<&str> = samples.iter().map(|s| s.as_str()).collect();
                        let tracks = loci::d4_tasks::prepare_tracks_from_file(
                            loci_args.infile.clone(),
                            Some(sample_refs),
                        )?;
                        let res = loci::d4_tasks::run_tasks_on_tracks(
                            tracks,
                            chrom_regions.clone(),
                            (loci_args.mean_depth_min, loci_args.mean_depth_max),
                            loci_args.depth_proportion,
                            loci_args.output_counts,
                        )?;
                        let fp = loci::write_d4::<PathBuf>(res, chroms.clone(), None)?;
                        temp_file_paths.push(fp);
                    }
                    loci::merge_d4_files(
                        loci_args.outfile.to_path_buf().into_std_path_buf(),
                        temp_file_paths,
                        population_map.get_popname_refs(),
                    )?;
                } else {
                    let tracks =
                        loci::d4_tasks::prepare_tracks_from_file(loci_args.infile.clone(), None)?;
                    let res = loci::d4_tasks::run_tasks_on_tracks(
                        tracks,
                        chrom_regions.clone(),
                        (loci_args.mean_depth_min, loci_args.mean_depth_max),
                        loci_args.depth_proportion,
                        loci_args.output_counts,
                    )?;
                    if bed_out {
                        loci::write_bed(loci_args.outfile, res, loci_args.output_counts)?;
                    } else {
                        loci::write_d4::<PathBuf>(
                            res,
                            chroms.clone(),
                            Some(loci_args.outfile.into_std_path_buf()),
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
