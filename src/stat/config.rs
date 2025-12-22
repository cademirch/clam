use crate::core::{
    contig::{validate_contig_consistency, Contig, ContigChunk, ContigSet},
    population::PopulationMap,
    zarr::CallableArrays,
};
use crate::stat::{
    roh::RohIndex,
    vcf::{build_vcf_reader, query::VcfQuery, variants::Ploidy},
    windows::WindowStrategy,
};
use color_eyre::eyre::WrapErr;
use color_eyre::Result;
use ndarray::Array2;
use noodles::core::{Position, Region};
use noodles::vcf::Header;
use std::{
    collections::HashSet,
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
    pub roh_bed_path: Option<PathBuf>,

    pub ploidy: Ploidy,
}

impl StatConfig {
    pub fn new(
        vcf_path: PathBuf,
        pop_file: Option<PathBuf>,
        callable_zarr: Option<PathBuf>,
        roh_bed: Option<PathBuf>,
        window_strategy: WindowStrategy,
        chunk_size: u64,
        include_contigs: Option<HashSet<String>>,
        exclude_contigs: Option<HashSet<String>>,
    ) -> Result<Self> {
        let (_, vcf_header) = build_vcf_reader(&vcf_path)?;

        let vcf_contigs = extract_contigs_from_vcf(&vcf_header)?;

        let vcf_samples: Vec<String> = vcf_header
            .sample_names()
            .iter()
            .map(|s| s.to_string())
            .collect();

        let pop_map = if let Some(pop_file_path) = pop_file {
            let pop_map = PopulationMap::from_file(pop_file_path)?;
            pop_map.validate_exact_match(&vcf_samples)?;
            pop_map
        } else {
            PopulationMap::default_from_samples(vcf_samples.clone())
        };

        let population_array = pop_map.membership_matrix(&vcf_samples)?;

        if let Some(ref zarr_path) = callable_zarr {
            let arrays = CallableArrays::open(zarr_path)?;
            let zarr_contigs = arrays.contigs();
            validate_contig_consistency(vec![
                (&vcf_path, vcf_contigs.clone()),
                (zarr_path, zarr_contigs.clone()),
            ])?;
        }
        //TODO: check zarr chunksize and use that
        let vcf_contigs = vcf_contigs.filter(include_contigs.as_ref(), exclude_contigs.as_ref());

        let ploidy = infer_ploidy(&vcf_path, &vcf_header)?;

        let chunks = vcf_contigs.to_chunks(chunk_size);

        Ok(Self {
            vcf_path,
            vcf_header,
            contig_set: vcf_contigs,
            chunks,
            pop_map,
            population_array,
            window_strategy,
            callable_zarr_path: callable_zarr,
            roh_bed_path: roh_bed,
            ploidy,
        })
    }

    /// Create a VcfQuery for processing a specific genomic chunk
    pub fn create_query(&self, chunk: &ContigChunk) -> Result<VcfQuery> {
        let start = Position::try_from((chunk.start+1) as usize).wrap_err("Invalid start position")?;
        let end = Position::try_from(chunk.end as usize).wrap_err("Invalid end position")?;

        let query_region = Region::new(chunk.contig_name.clone(), start..=end);

        let window_index = self.window_strategy.create_index(&query_region);

        let callable_loci = if let Some(ref zarr_path) = self.callable_zarr_path {
            Some(CallableArrays::open(zarr_path)?)
        } else {
            None
        };

        let roh_index = if let Some(ref bed_path) = self.roh_bed_path {
            Some(RohIndex::from_tabix_query(
                bed_path,
                &query_region,
                &self.vcf_header,
                self.pop_map.clone(),
            )?)
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
