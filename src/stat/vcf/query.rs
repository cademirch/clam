use crate::{
    core::{contig::ContigChunk, population::PopulationMap, zarr::CallableArrays},
    stat::{
        roh::RohIndex,
        vcf::variants::{genotypes, AlleleCounts, Ploidy},
        windows::{RegionExt, Window, WindowIdx, WindowIndex, WindowStats},
    },
};
use color_eyre::eyre::{bail, Result};
use ndarray::{Array1, Array2};
use noodles::core::Region;
use noodles::vcf::{self, Header, Record};
use rayon::prelude::*;
use std::path::PathBuf;

/// Helper to get variant position from record as usize since its nested in option<result<position>>
fn variant_start(record: &Record) -> Result<usize> {
    let pos = match record.variant_start() {
        Some(Ok(pos)) => pos,
        Some(Err(_)) | None => bail!(
            "Variant record: {:?} does not have a valid position",
            record
        ),
    };

    match pos.try_into() {
        Ok(u) => Ok(u),
        Err(_) => bail!(
            "Variant record: {:?} does not have a valid position",
            record
        ),
    }
}
pub struct VcfQuery {
    pub query_chunk: ContigChunk,
    pub query_region: Region,
    pub window_index: Box<dyn WindowIndex>,
    pub callable_loci: Option<CallableArrays>,
    pub roh_index: Option<RohIndex>,
    pub population_array: Array2<usize>,
    pub pop_map: PopulationMap,
    pub ploidy: Ploidy,
    pub vcf_path: PathBuf,
}
/// (window_idx, variant position, allele_counts)

pub type VariantResult = (usize, usize, AlleleCounts);
impl VcfQuery {
    /// Convert 1-based genomic position to 0-based array index within this query
    #[inline]
    fn genomic_to_array_idx(&self, position: usize) -> usize {
        position - self.query_region.start_usize()
    }
    fn process_variant_record(
        &self,
        record: &Record,
        header: &Header,
    ) -> Result<Option<VariantResult>> {
        let variant_position = variant_start(record)?;
        let pos_idx = variant_position - self.query_region.start_usize();
        let window_idx = match self.window_index.window_idx_for_pos(variant_position) {
            Some(idx) => idx,
            None => return Ok(None), // Variant outside all windows
        };
        let genotypes = genotypes(record, header)?;
        let samples_in_roh = self
            .roh_index
            .as_ref()
            .map(|roh_idx| roh_idx.samples_in_roh_at(variant_position));
        let allele_counts = AlleleCounts::from_genotypes(
            &genotypes,
            self.ploidy,
            samples_in_roh.as_deref(),
            &self.population_array,
        )?;

        Ok(Some((window_idx, variant_position, allele_counts)))
    }
    fn compute_invariant_refs(
        &self,
        position: usize, // 1-based genomic position
        callable_per_pop: &Array2<u16>,
    ) -> (Array1<usize>, Option<Array1<usize>>) {
        let array_idx = self.genomic_to_array_idx(position);

        let callable_at_pos = callable_per_pop.row(array_idx);
        let ref_alleles_per_sample = self.ploidy.as_usize();

        // All callable samples contribute ref alleles (invariant = all 0/0)
        let refs_all = callable_at_pos.mapv(|count| count as usize * ref_alleles_per_sample);

        // NON-ROH refs: callable - roh per population
        let refs_non_roh = self.roh_index.as_ref().map(|roh_idx| {
            let roh_counts_per_pop = roh_idx.roh_counts_per_pop(position);

            callable_at_pos
                .iter()
                .zip(roh_counts_per_pop.iter())
                .map(|(&callable, &roh)| {
                    callable.saturating_sub(roh) as usize * ref_alleles_per_sample
                })
                .collect()
        });

        (refs_all, refs_non_roh)
    }
    pub fn load_callable_data(&self) -> Result<Option<Array2<u16>>> {
        let res = match &self.callable_loci {
            Some(callable_arrays) => Some(
                callable_arrays
                    .read_chunk(&self.query_chunk.contig_name, self.query_chunk.chunk_idx)?,
            ),
            None => None,
        };
        Ok(res)
    }

    pub fn process(&self) -> Result<Vec<WindowStats>> {
        let (mut vcf_reader, vcf_header) = super::build_vcf_reader(&self.vcf_path)?;
        let variants = vcf_reader.query(&vcf_header, &self.query_region)?;

        let variant_stats: Vec<(WindowIdx, usize, AlleleCounts)> = variants
            .par_bridge()
            .filter_map(|result| {
                let record = result.ok()?;
                self.process_variant_record(&record, &vcf_header).ok()?
            })
            .collect();

        let mut windows: Vec<Window> = (0..self.window_index.num_windows())
            .map(|idx| {
                let region = self.window_index.window_region(idx);
                Window::new(region, self.pop_map.clone(), self.roh_index.is_some())
            })
            .collect();

        for (window_idx, position, counts) in variant_stats {
            windows[window_idx].add_variant(position, counts);
        }

        let callable_per_pop = self.load_callable_data()?;

        if let Some(ref callable) = callable_per_pop {
            let num_pops = self.pop_map.num_populations();

            windows.par_iter_mut().try_for_each(|window| {
                let window_start = window.region.start_usize();
                let window_end = window.region.end_usize();

                let mut total_refs = Array1::zeros(num_pops);
                let mut total_refs_non_roh = if self.roh_index.is_some() {
                    Some(Array1::zeros(num_pops))
                } else {
                    None
                };

                for position in window_start..=window_end {
                    if window.variant_positions.contains(&position) {
                        continue;
                    }

                    let (refs, refs_non_roh) = self.compute_invariant_refs(position, callable);

                    total_refs += &refs;
                    if let (Some(total_non_roh), Some(refs_non_roh)) =
                        (&mut total_refs_non_roh, refs_non_roh)
                    {
                        *total_non_roh += &refs_non_roh;
                    }
                }

                window.add_invariants(total_refs, total_refs_non_roh);

                Ok::<_, color_eyre::Report>(())
            })?;
        }

        let results = windows
            .into_par_iter()
            .map(|w| w.compute_stats())
            .collect::<Result<Vec<_>>>()?;

        Ok(results)
    }
}
