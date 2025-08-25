pub mod alleles;
pub mod callable;
pub mod chunks;
pub mod output;
pub mod windows;

use std::collections::HashSet;
use std::fs::File;
use std::num::NonZeroUsize;
use std::ops::Bound;
use std::path::Path;
use std::time;

// use alleles::process;
use anyhow::{anyhow, bail, Context, Result};
use bstr::{BString, ByteSlice};
use camino::Utf8PathBuf;
use clap::{ArgGroup, Parser};
use fnv::FnvHashMap;
use indexmap::IndexSet;
use log::{debug, info, trace, warn};
use noodles::bgzf::Reader;
use noodles::core::{Position, Region};
use noodles::csi::binning_index::ReferenceSequence;
use noodles::csi::BinningIndex;
use noodles::vcf::io::IndexedReader;
use noodles::vcf::variant::record::samples::keys::key;
use noodles::vcf::variant::record::samples::series::Value;
use noodles::vcf::variant::record::samples::Series;
use noodles::vcf::{self};
use windows::Window;

use crate::utils::{get_exclude_chromosomes, read_bed_regions, PopulationMapping};
use output::create_output_file;

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
    #[arg(long = "roh-file")]
    pub roh_file: Option<Utf8PathBuf>,

    /// Comma separated list of chromosomes to include (restrict analysis to)
    #[arg(short = 'i', value_delimiter = ',', num_args = 1.., conflicts_with("include_file"))]
    pub include: Option<Vec<String>>,

    /// Path to file with chromosomes to include, one per line
    #[arg(long = "include-file", conflicts_with("include"))]
    pub include_file: Option<Utf8PathBuf>,
}

pub fn build_vcf_reader(
    path: impl AsRef<Path>,
) -> Result<(IndexedReader<Reader<File>>, vcf::Header)> {
    let mut reader = vcf::io::indexed_reader::Builder::default()
        .build_from_path(path.as_ref())
        .context(format!(
            "Failed to read VCF file: {}",
            &path.as_ref().display()
        ))?;
    let header = reader
        .read_header()
        .with_context(|| "Failed to read VCF header")?;
    Ok((reader, header))
}

pub fn run_stat(args: StatArgs, progress_bar: Option<indicatif::ProgressBar>) -> Result<()> {
    let (header, ploidy) = get_vcf_header_and_ploidy(&args.vcf)?;

    let tbi_path = &args.vcf.with_extension("gz.tbi");
    let index = noodles::tabix::read(tbi_path)
        .context(format!("Failed to read tabix file: {}", tbi_path.as_str()))?;
    let index_header = index.header().context("Tabix file missing header")?;
    let index_seqs = index_header
        .reference_sequence_names()
        .into_iter()
        .map(|bs| String::from_utf8_lossy(&bs).into_owned())
        .collect();

    let pop_map = if let Some(pop_file) = &args.population_file {
        let file = File::open(&pop_file).with_context(|| {
            format!(
                "Failed to open population file at path: {}",
                pop_file.as_str()
            )
        })?;
        let pop_map = PopulationMapping::from_path(file, Some(header.sample_names()))?;
        pop_map.validate_sample_coverage(header.sample_names())?;
        pop_map
    } else {
        PopulationMapping::default(&header.sample_names())
    };

    let num_populations = pop_map.num_populations();
    let pop_names = pop_map.get_popname_refs();
    let exclude_chroms = get_exclude_chromosomes(&args.exclude, &args.exclude_file)?;
    let include_chroms = get_exclude_chromosomes(&args.include, &args.include_file)?;
    let seqlens = seqlens_vcf(
        &header,
        exclude_chroms.as_ref(),
        include_chroms.as_ref(),
        &index_seqs,
    )?;

    let regions = match (args.regions_file.clone(), args.window_size) {
        (Some(regions_file), _) => read_bed_regions(regions_file)?,
        (None, Some(window_size)) => regions_from_seqlens(window_size, seqlens, args.vcf.clone())?,
        _ => bail!("Either regions or windows are required!"),
    };

    let sites = match args.sites_file {
        Some(sites_file) => {
            let sites_regions = read_bed_regions(sites_file)?;
            let sites_map = sites_map(sites_regions)?;
            Some(make_region_sites_binary_search(&regions, sites_map)?)
        }
        None => None,
    };
    if let Some(callable_loci_file) = &args.callable_sites {
        if !callable_loci_file.exists() {
            bail!("Callable loci file {callable_loci_file} specified, but does not exist!")
        }
    }
    let windows = Window::from_regions(regions, &pop_map, sites, ploidy as u32);

    let worker_count = args.threads;

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

    let mut pi_writer = create_output_file(&outdir, "clam_pi.tsv")?;
    let mut dxy_writer = if num_populations > 1 {
        Some(create_output_file(&outdir, "clam_dxy.tsv")?)
    } else {
        None
    };
    let mut fst_writer = if num_populations > 1 {
        Some(create_output_file(&outdir, "clam_fst.tsv")?)
    } else {
        None
    };
    let mut het_writer = create_output_file(&outdir, "clam_het.tsv")?;

    for (_, window) in results.iter_mut().enumerate() {
        for pi_record in window.to_pi_records() {
            pi_writer.serialize(pi_record)?;
        }

        if let Some(ref mut dxy_writer) = dxy_writer {
            for dxy_record in window.to_dxy_records() {
                dxy_writer.serialize(dxy_record)?;
            }
        }

        if let Some(ref mut fst_writer) = fst_writer {
            for fst_record in window.to_fst_records() {
                fst_writer.serialize(fst_record)?;
            }
        }

        for het_record in window.to_het_records() {
            het_writer.serialize(het_record)?;
        }
    }

    pi_writer.flush()?;
    if let Some(ref mut dxy_writer) = dxy_writer {
        dxy_writer.flush()?;
    }
    het_writer.flush()?;
    Ok(())
}

pub fn regions_from_seqlens(
    window_size: usize,
    seqlens: FnvHashMap<String, usize>,
    vcf_path: Utf8PathBuf,
) -> Result<Vec<Region>> {
    let mut windows = vec![];
    let mut vcf_reader = vcf::io::indexed_reader::Builder::default().build_from_path(vcf_path)?;

    let header = vcf_reader.read_header()?;
    for (chrom, &length) in seqlens.iter() {
        // Calculate the number of windows for this chromosome based on its length
        let begin = Position::new(1).ok_or_else(|| anyhow!("Invalid start position"))?;
        let end = Position::new(length).ok_or_else(|| anyhow!("Invalid end position"))?;

        let region = Region::new(BString::from(chrom.as_str()), begin..=end);
        let time = time::Instant::now();
        let mut query = vcf_reader.query(&header, &region)?;

        let record = query.next();

        let first_variant_pos = if let Some(record) = record {
            let result = record?;
            result.variant_start().unwrap().unwrap().get()
        } else {
            1
        };
        // let first_variant_pos = query.count();

        // debug!(
        //     "Query count: {}. Time to query record count: {:?}",
        //     first_variant_pos,
        //     time.elapsed()
        // );
        let num_windows = ((length - first_variant_pos + 1) + window_size - 1) / window_size;

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

pub fn seqlens_vcf(
    header: &vcf::Header,
    exclude: Option<&HashSet<String>>,
    include: Option<&HashSet<String>>,
    index_seqs: &IndexSet<String>,
) -> Result<FnvHashMap<String, usize>> {
    let mut res = FnvHashMap::default();

    let contigs = header.contigs();
    if contigs.is_empty() {
        return Err(anyhow!(
            "No contig information found in the header. Please supply fasta index."
        ));
    }

    // Validate that all included contigs exist in the header
    if let Some(include_set) = include {
        for contig in include_set {
            if !contigs.contains_key(contig) {
                return Err(anyhow!(
                    "Included contig '{}' not found in VCF header",
                    contig
                ));
            }
        }
    }

    for (name, contig_map) in contigs {
        if let Some(length) = contig_map.length() {
            res.insert(name.clone(), length);
        } else {
            return Err(anyhow!("Contig {} has no length specified", name));
        }
    }

    res.retain(|key, _| index_seqs.contains(key));

    // Apply exclude filter if present
    if let Some(exclude_set) = exclude {
        if !exclude_set.is_empty() {
            res.retain(|contig, _| !exclude_set.contains(contig));
        }
    }

    // Apply include filter if present
    if let Some(include_set) = include {
        if !include_set.is_empty() {
            res.retain(|contig, _| include_set.contains(contig));
        }
    }

    // Check if any contigs remain after filtering
    if res.is_empty() {
        return Err(anyhow!(
            "No contigs remaining after applying include/exclude filters"
        ));
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
    let (mut reader, header) = build_vcf_reader(vcf_path.as_ref())?;

    let first_record = reader
        .records()
        .next()
        .transpose()
        .with_context(|| "Failed to read first VCF record")?
        .ok_or_else(|| anyhow::anyhow!("VCF file is empty or has no records"))?;

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
            Ok((None, _)) => ploidy += 1,            // Increment for missing alleles too
            Err(e) => return Err(anyhow!("Error parsing allele: {:?}", e)),
        }
    }

    warn!("Inferred ploidy: {} from VCF. If this is incorrect, specify ploidy with the option --ploidy.", ploidy);

    drop(reader);
    Ok((header.clone(), ploidy))
}

pub fn optimize_chunks<P: AsRef<Path>>(vcf_path: P, tbi_path: P) -> Result<()> {
    let (reader, header) = build_vcf_reader(vcf_path.as_ref())?;

    let index = noodles::tabix::read(tbi_path.as_ref()).context(format!(
        "Failed to read tabix file: {:?}",
        tbi_path.as_ref().display()
    ))?;
    let index_header = index.header().context("Tabix file missing header")?;
    let index_seqs = index_header.reference_sequence_names();

    for (idx, reference_seq) in index.reference_sequences().iter().enumerate() {
        if let Some(metadata) = reference_seq.metadata() {
            let mapped_records = metadata.mapped_record_count();
            let contig_name = index_seqs.get_index(idx).unwrap();
            debug!(
                "Contig: {}. Mapped record count: {}",
                contig_name, mapped_records
            );
        }
    }

    Ok(())
}
