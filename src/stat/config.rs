use crate::core::{
    contig::{validate_contig_consistency, Contig, ContigChunk, ContigSet},
    population::PopulationMap,
    utils::create_spinner,
    zarr::{CallableArrays, CallableLociType, SampleMaskArrays},
};
use crate::stat::{
    roh::RohIndex,
    vcf::{build_vcf_reader, query::VcfQuery, variants::Ploidy},
    windows::WindowStrategy,
};
use color_eyre::eyre::{bail, WrapErr};
use color_eyre::Result;
use log::warn;
use ndarray::Array2;
use noodles::core::{Position, Region};
use noodles::vcf::Header;
use std::{
    collections::HashSet,
    io,
    path::{Path, PathBuf},
};

pub struct StatConfig {
    pub vcf_path: PathBuf,
    pub vcf_header: Header,
    pub contig_set: ContigSet,
    pub chunks: Vec<ContigChunk>,
    pub pop_map: PopulationMap,
    pub population_array: Array2<usize>,
    pub window_strategy: WindowStrategy,

    // Store paths, not loaded data (load per-chunk)
    pub callable_zarr_path: Option<PathBuf>,
    pub callable_loci_type: Option<CallableLociType>,
    pub roh_bed_path: Option<PathBuf>,

    pub ploidy: Ploidy,

    /// Indices into VCF samples to filter genotypes (None = use all samples)
    pub sample_filter_indices: Option<Vec<usize>>,

    /// Samples included in the analysis (filtered if force_samples is enabled)
    pub analysis_samples: Vec<String>,
}

impl StatConfig {
    pub fn new(
        vcf_path: PathBuf,
        pop_map: Option<PopulationMap>,
        callable_zarr: Option<PathBuf>,
        roh_bed: Option<PathBuf>,
        window_strategy: WindowStrategy,
        chunk_size: u64,
        include_contigs: Option<HashSet<String>>,
        exclude_contigs: Option<HashSet<String>>,
        force_samples: bool,
    ) -> Result<Self> {
        let spinner = create_spinner("Opening VCF file...");
        let (_, vcf_header) = build_vcf_reader(&vcf_path)?;

        spinner.set_message("Extracting contigs from VCF header...");
        let vcf_contigs = extract_contigs_from_vcf(&vcf_header)?;

        let vcf_samples: Vec<String> = vcf_header
            .sample_names()
            .iter()
            .map(|s| s.to_string())
            .collect();

        let has_populations = pop_map.is_some();

        // Warn if force_samples is used without populations
        if force_samples && !has_populations {
            warn!("--force-samples has no effect without populations specified.");
        }

        spinner.set_message("Loading populations...");
        let pop_map = if let Some(pm) = pop_map {
            pm.validate_exact_match(&vcf_samples, force_samples)?;
            pm
        } else {
            PopulationMap::default_from_samples(vcf_samples.clone())
        };

        // Filter samples when force_samples is enabled with a population file
        let (analysis_samples, sample_filter_indices) = if force_samples && has_populations {
            let filtered = pop_map.filter_samples(&vcf_samples);
            if filtered.is_empty() {
                bail!(
                    "No samples from VCF found in population file. Cannot proceed with analysis."
                );
            }
            let indices = pop_map.filter_sample_indices(&vcf_samples);
            (filtered, Some(indices))
        } else {
            (vcf_samples, None)
        };

        let population_array = pop_map.membership_matrix(&analysis_samples)?;

        let callable_loci_type = if let Some(ref zarr_path) = callable_zarr {
            spinner.set_message("Validating callable sites zarr...");

            // Try to open as CallableArrays first to get metadata
            let arrays = CallableArrays::open(zarr_path)?;
            let zarr_contigs = arrays.contigs();
            validate_contig_consistency(vec![
                (&vcf_path, vcf_contigs.clone()),
                (zarr_path, zarr_contigs.clone()),
            ])?;

            // Detect callable loci type from metadata
            let loci_type = arrays
                .callable_loci_type()
                .unwrap_or(CallableLociType::PopulationCounts); // Backward compatibility

            // If it's SampleMasks, validate sample names match
            if loci_type == CallableLociType::SampleMasks {
                let mask_arrays = SampleMaskArrays::open(zarr_path)?;
                let zarr_samples = mask_arrays.column_names();

                // Verify all analysis samples are in the zarr
                for sample in &analysis_samples {
                    if !zarr_samples.contains(sample) {
                        bail!(
                            "Sample '{}' not found in callable loci zarr. \
                             Zarr contains: {:?}",
                            sample,
                            zarr_samples
                        );
                    }
                }
            }

            // If it's PopulationCounts, validate column count matches stat-time populations
            if loci_type == CallableLociType::PopulationCounts {
                let zarr_num_pops = arrays.column_names().len();
                let stat_num_pops = pop_map.num_populations();
                if zarr_num_pops != stat_num_pops {
                    bail!(
                        "Population count mismatch: callable zarr has {} population column(s) ({}) \
                         but stat is configured with {} population(s). \
                         Re-run 'clam loci' with matching populations, or omit -p/--samples \
                         to auto-read populations from the callable zarr.",
                        zarr_num_pops,
                        arrays.column_names().join(", "),
                        stat_num_pops,
                    );
                }
            }

            // Warn if using population counts with ROH data
            if loci_type == CallableLociType::PopulationCounts && roh_bed.is_some() {
                warn!(
                    "Using population-level callable counts with ROH data. \
                     Callable sites outside ROH will be approximated. \
                     For accurate per-sample heterozygosity, use per-sample callable masks."
                );
            }

            Some(loci_type)
        } else {
            None
        };

        if roh_bed.is_some() {
            spinner.set_message("Checking ROH BED file...");
        }

        //TODO: check zarr chunksize and use that
        let vcf_contigs = vcf_contigs.filter(include_contigs.as_ref(), exclude_contigs.as_ref());

        spinner.set_message("Inferring ploidy from VCF...");
        let ploidy = infer_ploidy(&vcf_path, &vcf_header)?;

        let chunks = vcf_contigs.to_chunks(chunk_size);

        spinner.finish_with_message(format!(
            "Loaded {} samples, {} contigs, {} chunks",
            analysis_samples.len(),
            vcf_contigs.len(),
            chunks.len()
        ));

        Ok(Self {
            vcf_path,
            vcf_header,
            contig_set: vcf_contigs,
            chunks,
            pop_map,
            population_array,
            window_strategy,
            callable_zarr_path: callable_zarr,
            callable_loci_type,
            roh_bed_path: roh_bed,
            ploidy,
            sample_filter_indices,
            analysis_samples,
        })
    }

    /// Create a VcfQuery for processing a specific genomic chunk
    pub fn create_query(&self, chunk: &ContigChunk) -> Result<VcfQuery> {
        let start =
            Position::try_from((chunk.start + 1) as usize).wrap_err("Invalid start position")?;
        let end = Position::try_from(chunk.end as usize).wrap_err("Invalid end position")?;

        let query_region = Region::new(chunk.contig_name.clone(), start..=end);

        let window_index = self.window_strategy.create_index(&query_region);

        let callable_loci = match (&self.callable_zarr_path, self.callable_loci_type) {
            (Some(zarr_path), Some(CallableLociType::SampleMasks)) => Some(
                super::vcf::query::CallableData::SampleMasks(SampleMaskArrays::open(zarr_path)?),
            ),
            (Some(zarr_path), _) => {
                // Default to PopulationCounts for backward compatibility
                Some(super::vcf::query::CallableData::PopulationCounts(
                    CallableArrays::open(zarr_path)?,
                ))
            }
            (None, _) => None,
        };

        let roh_index = if let Some(ref bed_path) = self.roh_bed_path {
            match RohIndex::from_tabix_query(
                bed_path,
                &query_region,
                &self.analysis_samples,
                self.pop_map.clone(),
            ) {
                Ok(index) => Some(index),
                Err(e) => {
                    // "missing reference sequence name" occurs when contig isn't in tabix index
                    if let Some(io_err) = e.root_cause().downcast_ref::<io::Error>() {
                        if io_err.kind() == io::ErrorKind::InvalidInput
                            && io_err
                                .to_string()
                                .contains("missing reference sequence name")
                        {
                            None
                        } else {
                            return Err(e);
                        }
                    } else {
                        return Err(e);
                    }
                }
            }
        } else {
            None
        };

        Ok(VcfQuery {
            query_chunk: chunk.clone(),
            query_region,
            window_index,
            callable_loci,
            roh_index,
            population_array: self.population_array.clone(),
            pop_map: self.pop_map.clone(),
            ploidy: self.ploidy,
            vcf_path: self.vcf_path.clone(),
            sample_filter_indices: self.sample_filter_indices.clone(),
            analysis_samples: self.analysis_samples.clone(),
        })
    }
}

fn extract_contigs_from_vcf(header: &Header) -> Result<ContigSet> {
    let contigs: Vec<Contig> = header
        .contigs()
        .iter()
        .map(|(name, map)| {
            let length = map.length().ok_or_else(|| {
                color_eyre::eyre::eyre!("Contig '{}' missing length in VCF header", name)
            })?;

            Ok(Contig::new(name.to_string(), length))
        })
        .collect::<Result<Vec<_>>>()?;

    Ok(ContigSet::new(contigs))
}

fn infer_ploidy(vcf_path: &Path, header: &Header) -> Result<Ploidy> {
    let (mut reader, _) = build_vcf_reader(vcf_path)?;

    // Read first record
    let mut record = noodles::vcf::variant::RecordBuf::default();
    reader
        .read_record_buf(header, &mut record)
        .wrap_err("Failed to read first VCF record")?;

    Ploidy::from_record_buf(&record, header).wrap_err("Failed to infer ploidy from first variant")
}
