use crate::{
    core::{
        contig::ContigChunk,
        population::PopulationMap,
        zarr::{CallableArrays, SampleMaskArrays},
    },
    stat::{
        roh::RohIndex,
        vcf::variants::{genotypes, AlleleCounts, Ploidy},
        windows::{RegionExt, Window, WindowIdx, WindowIndex, WindowStats},
    },
};
use color_eyre::eyre::{bail, Result};
use log::debug;
use ndarray::{Array1, Array2};
use noodles::core::Region;
use noodles::vcf::{Header, Record};
use rayon::prelude::*;
use std::{io, path::PathBuf};

/// Callable loci data - either per-population counts or per-sample masks
pub enum CallableData {
    PopulationCounts(CallableArrays),
    SampleMasks(SampleMaskArrays),
}

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
    pub callable_loci: Option<CallableData>,
    pub roh_index: Option<RohIndex>,
    pub population_array: Array2<usize>,
    pub pop_map: PopulationMap,
    pub ploidy: Ploidy,
    pub vcf_path: PathBuf,
    /// Indices into VCF samples to filter genotypes (None = use all samples)
    pub sample_filter_indices: Option<Vec<usize>>,
    /// Sample names for the analysis (needed for sample mask lookups)
    pub analysis_samples: Vec<String>,
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
        let window_idx = match self.window_index.window_idx_for_pos(variant_position) {
            Some(idx) => idx,
            None => return Ok(None), // Variant outside all windows
        };
        let mut gts = genotypes(record, header)?;

        // Filter genotypes to only include analysis samples
        if let Some(ref indices) = self.sample_filter_indices {
            gts = indices.iter().map(|&i| gts[i]).collect();
        }

        let samples_in_roh = self
            .roh_index
            .as_ref()
            .map(|roh_idx| roh_idx.samples_in_roh_at(variant_position));
        let allele_counts = AlleleCounts::from_genotypes(
            &gts,
            self.ploidy,
            samples_in_roh.as_deref(),
            &self.population_array,
        )?;

        Ok(Some((window_idx, variant_position, allele_counts)))
    }

    /// Load callable data as population counts (u16 per population)
    fn load_callable_pop_counts(&self) -> Result<Option<Array2<u16>>> {
        match &self.callable_loci {
            Some(CallableData::PopulationCounts(arrays)) => {
                let data = arrays.read_chunk(
                    &self.query_chunk.contig_name,
                    self.query_chunk.chunk_idx,
                    None,
                )?;
                Ok(Some(data))
            }
            _ => Ok(None),
        }
    }

    /// Load callable data as per-sample masks (bool per sample)
    fn load_callable_sample_masks(&self) -> Result<Option<Array2<bool>>> {
        match &self.callable_loci {
            Some(CallableData::SampleMasks(arrays)) => {
                let data = arrays.read_chunk(
                    &self.query_chunk.contig_name,
                    self.query_chunk.chunk_idx,
                    None,
                )?;
                Ok(Some(data))
            }
            _ => Ok(None),
        }
    }

    /// Build a mapping from analysis sample indices to zarr column indices
    fn build_sample_to_zarr_idx(&self) -> Option<Vec<usize>> {
        match &self.callable_loci {
            Some(CallableData::SampleMasks(arrays)) => {
                let zarr_samples = arrays.column_names();
                let mapping: Vec<usize> = self
                    .analysis_samples
                    .iter()
                    .map(|sample| {
                        zarr_samples
                            .iter()
                            .position(|s| s == sample)
                            .expect("Sample should exist in zarr (validated at config time)")
                    })
                    .collect();
                Some(mapping)
            }
            _ => None,
        }
    }

    pub fn process(&self) -> Result<Vec<WindowStats>> {
        let (mut vcf_reader, vcf_header) = super::build_vcf_reader(&self.vcf_path)?;
        let variants = match vcf_reader.query(&vcf_header, &self.query_region) {
            Ok(v) => v,
            Err(e) => {
                // Handle case where contig exists in VCF header but not in tabix index
                // This can happen when the VCF has no variants on that contig
                if e.kind() == io::ErrorKind::InvalidInput
                    && e.to_string()
                        .contains("region reference sequence does not exist")
                {
                    debug!("Skipping region {} (not in VCF index)", self.query_region);
                    return Ok(vec![]);
                }
                return Err(e.into());
            }
        };

        let variant_stats: Vec<(WindowIdx, usize, AlleleCounts)> = variants
            .par_bridge()
            .filter_map(|result| {
                let record = result.ok()?;
                self.process_variant_record(&record, &vcf_header).ok()?
            })
            .collect();

        let num_samples = self.analysis_samples.len();
        let has_sample_masks = matches!(self.callable_loci, Some(CallableData::SampleMasks(_)));

        let mut windows: Vec<Window> = (0..self.window_index.num_windows())
            .map(|idx| {
                let region = self.window_index.window_region(idx);
                Window::new(
                    region,
                    self.pop_map.clone(),
                    self.roh_index.is_some(),
                    num_samples,
                    has_sample_masks,
                )
            })
            .collect();

        for (window_idx, position, counts) in variant_stats {
            windows[window_idx].add_variant(position, counts);
        }

        // Process callable sites based on the type of callable data available
        match &self.callable_loci {
            Some(CallableData::SampleMasks(_)) => {
                self.process_callable_sample_masks(&mut windows)?;
            }
            Some(CallableData::PopulationCounts(_)) => {
                self.process_callable_pop_counts(&mut windows)?;
            }
            None => {
                // No callable data - skip invariant site processing
            }
        }

        let results = windows
            .into_par_iter()
            .map(|w| w.compute_stats(&self.analysis_samples))
            .collect::<Result<Vec<_>>>()?;

        Ok(results)
    }

    /// Process callable sites using per-sample masks
    fn process_callable_sample_masks(&self, windows: &mut [Window]) -> Result<()> {
        let callable_masks = self.load_callable_sample_masks()?.unwrap();
        let sample_to_zarr = self.build_sample_to_zarr_idx().unwrap();
        let num_samples = self.analysis_samples.len();
        let num_pops = self.pop_map.num_populations();
        let num_pop_pairs = if num_pops >= 2 {
            num_pops * (num_pops - 1) / 2
        } else {
            0
        };

        windows.par_iter_mut().try_for_each(|window| {
            let window_start = window.region.start_usize();
            let window_end = window.region.end_usize();

            let mut callable_per_sample = Array1::<usize>::zeros(num_samples);
            let mut callable_per_sample_not_in_roh = Array1::<usize>::zeros(num_samples);

            // For pi/dxy calculations
            let mut total_comparisons = Array1::zeros(num_pops);
            let mut total_dxy_comparisons = vec![0; num_pop_pairs];

            for position in window_start..=window_end {
                if window.variant_positions.contains(&position) {
                    continue;
                }

                let array_idx = self.genomic_to_array_idx(position);
                let masks_at_pos = callable_masks.row(array_idx);

                // Get samples in ROH at this position
                let samples_in_roh: Vec<usize> = self
                    .roh_index
                    .as_ref()
                    .map(|idx| idx.samples_in_roh_at(position))
                    .unwrap_or_default();

                // Track per-sample callable sites
                for (sample_idx, &zarr_idx) in sample_to_zarr.iter().enumerate() {
                    if masks_at_pos[zarr_idx] {
                        callable_per_sample[sample_idx] += 1;

                        if !samples_in_roh.contains(&sample_idx) {
                            callable_per_sample_not_in_roh[sample_idx] += 1;
                        }
                    }
                }

                // Count callable samples per population for pi/dxy
                let mut callable_per_pop = vec![0usize; num_pops];

                for (sample_idx, &zarr_idx) in sample_to_zarr.iter().enumerate() {
                    if masks_at_pos[zarr_idx] {
                        if let Some(sample_name) = self.analysis_samples.get(sample_idx) {
                            if let Some((pop_idx, _)) = self.pop_map.lookup(sample_name) {
                                callable_per_pop[pop_idx] += 1;
                            }
                        }
                    }
                }

                // Compute comparisons for pi
                for pop_idx in 0..num_pops {
                    let haps = callable_per_pop[pop_idx] * self.ploidy.as_usize();
                    let comps = if haps > 1 { (haps * (haps - 1)) / 2 } else { 0 };
                    total_comparisons[pop_idx] += comps;
                }

                // Compute comparisons for dxy
                if num_pops >= 2 {
                    let mut pair_idx = 0;
                    for i in 0..num_pops {
                        let haps_i = callable_per_pop[i] * self.ploidy.as_usize();
                        for j in (i + 1)..num_pops {
                            let haps_j = callable_per_pop[j] * self.ploidy.as_usize();
                            total_dxy_comparisons[pair_idx] += haps_i * haps_j;
                            pair_idx += 1;
                        }
                    }
                }
            }

            window.add_callable_sites(callable_per_sample, Some(callable_per_sample_not_in_roh));

            window.add_invariant_comparisons(total_comparisons, total_dxy_comparisons);

            Ok::<_, color_eyre::Report>(())
        })?;

        Ok(())
    }

    /// Process callable sites using per-population counts
    fn process_callable_pop_counts(&self, windows: &mut [Window]) -> Result<()> {
        let callable_per_pop = self.load_callable_pop_counts()?.unwrap();
        let num_pops = self.pop_map.num_populations();
        let num_pop_pairs = if num_pops >= 2 {
            num_pops * (num_pops - 1) / 2
        } else {
            0
        };

        windows.par_iter_mut().try_for_each(|window| {
            let window_start = window.region.start_usize();
            let window_end = window.region.end_usize();

            let mut total_comparisons = Array1::zeros(num_pops);
            let mut total_dxy_comparisons = vec![0; num_pop_pairs];

            // For population-level heterozygosity (summed approach)
            let mut callable_per_pop_total = Array1::<usize>::zeros(num_pops);
            let mut callable_per_pop_not_in_roh = if self.roh_index.is_some() {
                Some(Array1::<usize>::zeros(num_pops))
            } else {
                None
            };

            for position in window_start..=window_end {
                if window.variant_positions.contains(&position) {
                    continue;
                }

                let array_idx = self.genomic_to_array_idx(position);
                let callable_at_pos = callable_per_pop.row(array_idx);

                // Get ROH counts per population if available
                let roh_counts: Vec<u16> = self
                    .roh_index
                    .as_ref()
                    .map(|idx| idx.roh_counts_per_pop(position))
                    .unwrap_or_else(|| vec![0; num_pops]);

                for pop_idx in 0..num_pops {
                    let callable_samples = callable_at_pos[pop_idx] as usize;
                    callable_per_pop_total[pop_idx] += callable_samples;

                    let haps = callable_samples * self.ploidy.as_usize();
                    let comps = if haps > 1 { (haps * (haps - 1)) / 2 } else { 0 };
                    total_comparisons[pop_idx] += comps;

                    // Approximate non-ROH callable (population counts mode)
                    if let Some(ref mut non_roh_total) = callable_per_pop_not_in_roh {
                        let non_roh_samples = (callable_at_pos[pop_idx])
                            .saturating_sub(roh_counts[pop_idx])
                            as usize;
                        non_roh_total[pop_idx] += non_roh_samples;
                    }
                }

                if num_pops >= 2 {
                    let mut pair_idx = 0;
                    for i in 0..num_pops {
                        let haps_i = callable_at_pos[i] as usize * self.ploidy.as_usize();
                        for j in (i + 1)..num_pops {
                            let haps_j = callable_at_pos[j] as usize * self.ploidy.as_usize();
                            total_dxy_comparisons[pair_idx] += haps_i * haps_j;
                            pair_idx += 1;
                        }
                    }
                }
            }

            window.add_callable_sites_per_pop(callable_per_pop_total, callable_per_pop_not_in_roh);

            window.add_invariant_comparisons(total_comparisons, total_dxy_comparisons);

            Ok::<_, color_eyre::Report>(())
        })?;

        Ok(())
    }
}
