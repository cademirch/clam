pub mod alleles;
pub mod callable;
pub mod windows;

use crate::utils::{get_exclude_chromosomes, read_bed_regions, PopulationMapping};
// use alleles::process;
use anyhow::{anyhow, bail, Context, Result};
use bstr::{BString, ByteSlice};
use camino::Utf8PathBuf;
use clap::{ArgGroup, Parser};
use fnv::FnvHashMap;
use log::{info, trace, warn};
use noodles::core::{Position, Region};
use noodles::csi::BinningIndex;
use noodles::vcf::{
    self,
    variant::record::samples::{keys::key, series::Value, Series},
};
use std::collections::HashSet;
use std::fs::File;
use std::num::NonZeroUsize;
use std::ops::Bound;
use std::path::Path;
use windows::Window;


#[derive(Parser, Debug, Clone)]
#[command(about = "Calculate population genetic statistics from VCF using callable sites.")]
#[command(group(ArgGroup::new("regions")
    .required(true)
    .multiple(false)
    .args(&["window_size", "regions_file"])),
    group(ArgGroup::new("exclude chromosomes")
    .required(false)
    .multiple(false)
    .args(&["exclude", "exclude_file"])))]
pub struct StatArgs {
    /// Path to input VCF file
    #[arg(required = true)]
    pub vcf: Utf8PathBuf,

    /// Path to input callable sites D4 file from clam loci
    pub callable_sites: Option<Utf8PathBuf>,

    /// Where to write output files. Defaults to current working directory.
    #[arg(short = 'o', long = "outdir")]
    pub outdir: Option<Utf8PathBuf>,

    /// Number of threads to use
    #[arg(short = 't', long = "threads", default_value_t = NonZeroUsize::new(1).unwrap())]
    pub threads: NonZeroUsize,

    /// Size of windows for statistics in bp. Conflicts with 'regions-file'
    #[arg(short = 'w', long = "window-size")]
    pub window_size: Option<usize>,

    /// File specifying regions to calculate statistics for. Conflicts with 'window-size'
    #[arg(short = 'r', long = "regions-file")]
    pub regions_file: Option<Utf8PathBuf>,

    /// Specify sites to consider for calculations. Bed format.
    #[arg(short = 's', long = "sites-file")]
    pub sites_file: Option<Utf8PathBuf>,

    /// Path to file that defines populations. Tab separated: sample, population_name
    #[arg(short = 'p', long = "population-file")]
    pub population_file: Option<Utf8PathBuf>,

    /// Path to fasta index for reference VCF was called against. Only needed if VCF does not have contig info in the header.
    #[arg(short = 'f', long = "fai")]
    pub fasta_index: Option<Utf8PathBuf>,

    /// Comma separated list of chromosomes to exclude
    #[arg(short = 'x', value_delimiter = ',', num_args = 1..)]
    pub exclude: Option<Vec<String>>,

    /// Path to file with chromosomes to exclude, one per line
    #[arg(long = "exclude-file", conflicts_with("exclude"))]
    pub exclude_file: Option<Utf8PathBuf>,

    /// Path to RoH file.
    #[arg(long = "roh-file", conflicts_with("callable_sites"))]
    pub roh_file: Option<Utf8PathBuf>,
}

#[derive(serde::Serialize)]
struct PiRecord {
    population_name: String,
    chrom: String,
    start: u32,
    end: u32,
    pi: f32,
    comparisons: u32,
    differences: u32,
}
#[derive(serde::Serialize)]
struct DxyRecord {
    population1_name: String,
    population2_name: String,
    chrom: String,
    start: u32,
    end: u32,
    dxy: f32,
    comparisons: u32,
    differences: u32,
}

pub fn run_stat(args: StatArgs, progress_bar: Option<indicatif::ProgressBar>) -> Result<()> {
    // Load VCF header and determine ploidy
    let (header, ploidy) = get_vcf_header_and_ploidy(&args.vcf)?;
    
    let tbi_path = &args.vcf.with_extension("gz.tbi");
    if !tbi_path.exists() {
        bail!("Couldn't find tabix index: {}", tbi_path);
    }
    let index = noodles::tabix::read(tbi_path)?;
    let index_header = index.header().context("Tabix file missing header")?;
    let index_seqs = index_header.reference_sequence_names();

    let pop_map = if let Some(pop_file) = &args.population_file {
        PopulationMapping::from_path(pop_file, Some(header.sample_names()))?
    } else {
        // Default mapping: 1 population, all samples in the same population
        PopulationMapping::default(&header)
    };
    
    // Set up sample-to-population mapping, if population file is provided
    let num_populations = pop_map.num_populations;
    
    let pop_names = pop_map.get_popname_refs();

    // Get sequence lengths for regions based on VCF header
    let mut seqlens = seqlens_vcf(&header)?;
    seqlens.retain(|key, _| index_seqs.contains(key));
    let exclude_chroms = get_exclude_chromosomes(&args.exclude, &args.exclude_file)?;

    // Modify `seqlens` in place if `exclude_chroms` has elements
    if !exclude_chroms.is_empty() {
        seqlens.retain(|contig, _| !exclude_chroms.contains(contig));
    }

    let regions = if let Some(regions_file) = args.regions_file.clone() {
        read_bed_regions(regions_file)?
    } else if let Some(window_size) = args.window_size {
        regions_from_seqlens(window_size, seqlens)?
    } else {
        bail!("Either regions or windows are required!")
    };

    let sites = if let Some(sites_file) = args.sites_file {
        let sites_regions = read_bed_regions(sites_file)?;
        let sites_mapp = sites_map(sites_regions)?;
        Some(make_region_sites_binary_search(&regions, sites_mapp)?)
    } else {
        None
    };

    // Set up windows from regions
    let windows = Window::from_regions(regions, &pop_map, sites, ploidy as u32);

    // Define worker count from threads argument
    let worker_count = args.threads;

    // Process the VCF data
    info!("Starting VCF processing...");
    let mut results = windows::process_windows(
        args.vcf,
        args.callable_sites,
        args.roh_file,
        worker_count,
        windows,
        progress_bar,
        pop_map.clone(),
    )?;

    let outdir = args.outdir.unwrap_or_else(|| Utf8PathBuf::from("."));

    // Create the clam_pi.tsv file
    let clam_pi_path = outdir.join("clam_pi.tsv");
    let clam_pi_file = File::create(&clam_pi_path)?;
    let mut clam_pi_writer = csv::WriterBuilder::new()
        .delimiter(b'\t')
        .from_writer(clam_pi_file);
    // Conditionally create the clam_dxy.tsv file
    let mut clam_dxy_writer = if num_populations >= 2 {
        let clam_dxy_path = outdir.join("clam_dxy.tsv");
        let clam_dxy_file = File::create(&clam_dxy_path)?;
        let clam_dxy_writer = csv::WriterBuilder::new()
            .delimiter(b'\t')
            .from_writer(clam_dxy_file);
        Some(clam_dxy_writer)
    } else {
        None
    };

    // single threaded version for both
    // threaded too

    for window in results.iter_mut() {
        for (pop_idx, population) in window.populations.iter().enumerate() {
            let comps = population.within_comps;
            let diffs = population.within_diffs;
            let pi = if comps > 0 {
                diffs as f32 / comps as f32
            } else {
                0.0 // Handle zero comparisons
            };
            let pop_names = if let Some(owned_pop_names) = pop_names.clone() {
                // owned_pop_names is now a fully owned Vec<&str>, cloned from outer pop_names
                owned_pop_names
            } else {
                vec!["pop1"]
            };
            let (chrom, begin, end) = window.get_region_info();
            clam_pi_writer.serialize(PiRecord {
                population_name: pop_names[pop_idx].to_string(),
                chrom: chrom.to_string(),
                start: begin,
                end: end,
                pi: pi,
                comparisons: comps,
                differences: diffs,
            })?;
        }
        if let Some(dxy_writer) = clam_dxy_writer.as_mut() {
            let pop_names = if let Some(owned_pop_names) = pop_names.clone() {
                // owned_pop_names is now a fully owned Vec<&str>, cloned from outer pop_names
                owned_pop_names
            } else {
                vec!["pop1"]
            };
            let (chrom, begin, end) = window.get_region_info();
            for pop1 in 0..num_populations {
                for pop2 in (pop1 + 1)..num_populations {
                    // Get the index for this population pair
                    let pair_idx = window.get_pair_index(pop1, pop2);

                    // Retrieve comparisons and differences for DXY calculation
                    if let (Some(comps), Some(diffs)) = (
                        window.dxy_comps.as_ref().and_then(|c| c.get(pair_idx)),
                        window.dxy_diffs.as_ref().and_then(|d| d.get(pair_idx)),
                    ) {
                        let dxy_value = if *comps > 0 {
                            *diffs as f32 / *comps as f32
                        } else {
                            0.0 // Handle zero comparisons
                        };

                        dxy_writer.serialize(DxyRecord {
                            population1_name: pop_names[pop1].to_string(),
                            population2_name: pop_names[pop2].to_string(),
                            chrom: chrom.to_string(),
                            start: begin,
                            end: end,
                            dxy: dxy_value,
                            comparisons: *comps,
                            differences: *diffs,
                        })?;
                    }
                }
            }
        }
    }
    Ok(())
}

pub fn regions_from_seqlens(
    window_size: usize,
    seqlens: FnvHashMap<String, usize>,
) -> Result<Vec<Region>> {
    let mut windows = vec![];
    for (chrom, &length) in seqlens.iter() {
        // Calculate the number of windows for this chromosome based on its length
        let num_windows = (length + window_size - 1) / window_size;

        for window_idx in 0..num_windows {
            // Calculate the 1-based start and end positions for each window
            let begin = Position::new(window_idx * window_size + 1)
                .ok_or_else(|| anyhow!("Invalid start position"))?;
            let end = Position::new(((window_idx + 1) * window_size).min(length))
                .ok_or_else(|| anyhow!("Invalid end position"))?;

            // Create a new region with the chromosome name and bounds
            let region = Region::new(BString::from(chrom.as_str()), begin..=end);
            windows.push(region);
        }
    }

    Ok(windows)
}

pub fn seqlens_vcf(header: &vcf::Header) -> Result<FnvHashMap<String, usize>> {
    let mut res = FnvHashMap::default();

    let contigs = header.contigs();
    if contigs.is_empty() {
        return Err(anyhow!(
            "No contig information found in the header. Please supply fasta index."
        ));
    }

    for (name, contig_map) in contigs {
        if let Some(length) = contig_map.length() {
            res.insert(name.clone(), length);
        } else {
            return Err(anyhow!("Contig {} has no length specified", name));
        }
    }

    Ok(res)
}
fn make_region_sites_binary_search(
    regions: &[Region],
    mut sites: FnvHashMap<String, Vec<u32>>,
) -> Result<Vec<HashSet<u32>>> {
    let mut region_sites = Vec::with_capacity(regions.len());

    // Ensure all chromosome sites are sorted
    for chrom_sites in sites.values_mut() {
        chrom_sites.sort_unstable(); // Sort each list of positions in place
    }

    for region in regions {
        let mut site_set = HashSet::new();

        if let Some(chrom_sites) = sites.get(
            region
                .name()
                .to_str()
                .map_err(|_| anyhow!("Failed to convert region name to string"))?,
        ) {
            let (start, end) = get_region_positions(region)?;
            let start_idx = chrom_sites.partition_point(|&pos| pos < start as u32);
            let end_idx = chrom_sites.partition_point(|&pos| pos < end as u32);

            // Collect positions within the range
            site_set.extend(&chrom_sites[start_idx..end_idx]);
        }

        region_sites.push(site_set);
    }

    Ok(region_sites)
}

/// Extracts the start and end positions from a `Region`, normalizing bounds to inclusive positions.
fn get_region_positions(region: &Region) -> Result<(usize, usize)> {
    let start = match region.start() {
        Bound::Included(pos) => pos.get(),
        Bound::Excluded(pos) => pos.get() + 1,
        Bound::Unbounded => anyhow::bail!("Start bound of region is unbounded"),
    };

    let end = match region.end() {
        Bound::Included(pos) => pos.get(),
        Bound::Excluded(pos) => pos.get() - 1,
        Bound::Unbounded => anyhow::bail!("End bound of region is unbounded"),
    };

    Ok((start, end))
}

pub fn sites_map(sites_regions: Vec<Region>) -> Result<FnvHashMap<String, Vec<u32>>> {
    let mut res: FnvHashMap<String, Vec<u32>> = FnvHashMap::default();

    for region in sites_regions {
        let name = String::from_utf8(region.name().to_vec())
            .map_err(|e| anyhow::anyhow!("Invalid UTF-8 in region name: {}", e))?;

        let (start, end) = get_region_positions(&region)?;

        let positions = res.entry(name).or_insert_with(Vec::new);
        for pos in start..=end {
            positions.push(pos as u32);
        }
    }

    Ok(res)
}
pub fn get_vcf_header_and_ploidy<P: AsRef<Path>>(vcf_path: P) -> Result<(vcf::Header, usize)> {
    // Create a VCF reader from the provided path
    let mut reader = vcf::io::reader::Builder::default()
        .build_from_path(vcf_path.as_ref())
        .with_context(|| format!("Failed to open VCF file: {:?}", vcf_path.as_ref()))?;

    // Read the header
    let header = reader
        .read_header()
        .with_context(|| "Failed to read VCF header")?;

    // Read the first record to determine ploidy
    let first_record = reader
        .records()
        .next()
        .transpose()
        .with_context(|| "Failed to read first VCF record")?
        .ok_or_else(|| anyhow::anyhow!("VCF file is empty or has no records"))?;

    // Extract the genotype (GT) field and count alleles
    let samples = first_record.samples();
    let gt_series = samples
        .select(key::GENOTYPE)
        .ok_or_else(|| anyhow!("Malformed variant record: {:?}", first_record))?;

    let Some(Value::Genotype(genotype)) = gt_series
        .iter(&header)
        .next()
        .context("Failed to get genotype.")??
    else {
        return Err(anyhow!(
            "GT field is missing or invalid in the first record"
        ));
    };

    let mut ploidy = 0;
    for result in genotype.iter() {
        match result {
            Ok((Some(_position), _)) => ploidy += 1, // Increment for each non-missing allele
            Ok((None, _)) => ploidy += 1,               // Increment for missing alleles too
            Err(e) => return Err(anyhow!("Error parsing allele: {:?}", e)),
        }
    }
    
    warn!("Inferred ploidy: {} from VCF. If this is incorrect, specify ploidy with the option --ploidy.", ploidy);
    
    drop(reader);
    Ok((header.clone(), ploidy))
}

