use crate::core::population::PopulationMap;
use crate::stat::vcf::variants::AlleleCounts;
use bstr::ByteSlice;
use color_eyre::eyre::Ok;
use color_eyre::Result;
use ndarray::Array1;
use noodles::core::{Position, Region};
use rust_lapper::{Interval, Lapper};
use serde::Serialize;
use std::ops::Bound;

use std::collections::HashSet;

#[derive(Debug, Clone, Serialize)]
pub struct SerializableRegion {
    pub chrom: String,
    pub start: usize,
    pub end: usize,
}

impl From<&Region> for SerializableRegion {
    fn from(region: &Region) -> Self {
        Self {
            chrom: region.name_as_str().to_string(),
            start: region.start_usize(),
            end: region.end_usize(),
        }
    }
}
pub trait RegionExt {
    fn start_usize(&self) -> usize;
    fn end_usize(&self) -> usize;
    fn name_as_str(&self) -> &str;
}
impl RegionExt for Region {
    fn start_usize(&self) -> usize {
        match self.start() {
            Bound::Included(pos) => pos.get().try_into().unwrap(),
            Bound::Excluded(_) => panic!("Region starts should always be Included."),
            Bound::Unbounded => panic!("Region starts should always be bounded."),
        }
    }

    fn end_usize(&self) -> usize {
        match self.end() {
            Bound::Included(pos) => pos.get().try_into().unwrap(),
            Bound::Excluded(_) => panic!("Region ends should always be Included."),
            Bound::Unbounded => panic!("Region ends should always be bounded."),
        }
    }

    fn name_as_str(&self) -> &str {
        self.name()
            .to_str()
            .expect("Region name is not valid UTF-8")
    }
}
pub type WindowIdx = usize;
pub trait WindowIndex: Send + Sync {
    /// Get the window index for a position
    fn window_idx_for_pos(&self, pos: usize) -> Option<WindowIdx>;

    /// Get the region for a specific window index
    fn window_region(&self, idx: usize) -> Region;

    /// Number of windows
    fn num_windows(&self) -> usize;
}

/// O(1) lookup for fixed-size windows - computes regions on demand
pub struct FixedWindowIndex {
    chrom: String,
    region_start: usize,
    region_end: usize,
    window_size: usize,
    num_windows: usize,
}

impl FixedWindowIndex {
    pub fn new(region: &Region, window_size: usize) -> Self {
        let start = region.start_usize();
        let end = region.end_usize();
        let chrom = region.name_as_str().to_string();

        let num_windows = (end - start).div_ceil(window_size);

        Self {
            chrom,
            region_start: start,
            region_end: end,
            window_size,
            num_windows,
        }
    }
}

impl WindowIndex for FixedWindowIndex {
    #[inline]
    fn window_idx_for_pos(&self, pos: usize) -> Option<usize> {
        if pos < self.region_start || pos >= self.region_end {
            return None;
        }
        let idx = (pos - self.region_start) / self.window_size;
        if idx < self.num_windows {
            Some(idx)
        } else {
            None
        }
    }

    fn window_region(&self, idx: usize) -> Region {
        let start = self.region_start + (idx * self.window_size);
        // Use 1-based closed coordinates: window 0 is [1, 10000], window 1 is [10001, 20000], etc.
        // This matches the old clam behavior and standard genomic conventions.
        let end = (start + self.window_size - 1).min(self.region_end);
        let start_pos = Position::try_from(start).expect("Valid position");
        let end_pos = Position::try_from(end).expect("Valid position");

        Region::new(&*self.chrom, start_pos..=end_pos)
    }

    fn num_windows(&self) -> usize {
        self.num_windows
    }
}

pub struct RegionWindowIndex {
    windows: Vec<Region>,
    lapper: Lapper<usize, usize>,
}

impl RegionWindowIndex {
    pub fn new(region: &Region, all_windows: &[Region]) -> Self {
        let windows: Vec<Region> = all_windows
            .iter()
            .filter(|w| {
                w.name_as_str() == region.name_as_str()
                    && w.start_usize() < region.end_usize()
                    && w.end_usize() > region.start_usize()
            })
            .cloned()
            .collect();

        let intervals: Vec<Interval<usize, usize>> = windows
            .iter()
            .enumerate()
            .map(|(idx, w)| Interval {
                start: w.start_usize(),
                stop: w.end_usize(),
                val: idx,
            })
            .collect();

        Self {
            windows,
            lapper: Lapper::new(intervals),
        }
    }
}

impl WindowIndex for RegionWindowIndex {
    #[inline]
    fn window_idx_for_pos(&self, pos: usize) -> Option<usize> {
        self.lapper.find(pos, pos + 1).next().map(|iv| iv.val)
    }

    fn window_region(&self, idx: usize) -> Region {
        self.windows[idx].clone()
    }

    fn num_windows(&self) -> usize {
        self.windows.len()
    }
}

/// Strategy for creating windows
pub enum WindowStrategy {
    FixedSize(usize),
    Regions(Vec<Region>), // bed supplied windows
}

impl WindowStrategy {
    pub fn create_index(&self, region: &Region) -> Box<dyn WindowIndex> {
        match self {
            WindowStrategy::FixedSize(size) => Box::new(FixedWindowIndex::new(region, *size)),
            WindowStrategy::Regions(windows) => Box::new(RegionWindowIndex::new(region, windows)),
        }
    }
}

pub struct Window {
    pub region: Region,
    pop_map: PopulationMap,
    pub variant_positions: HashSet<usize>,

    // Per-sample heterozygosity tracking (for all samples)
    pub het_counts: Array1<usize>,
    pub het_counts_not_in_roh: Option<Array1<usize>>,

    // Per-sample callable site tracking (only when sample masks available)
    pub callable_sites: Option<Array1<usize>>,
    pub callable_sites_not_in_roh: Option<Array1<usize>>,

    // Per-population callable site tracking (for population counts mode)
    pub callable_sites_per_pop: Option<Array1<usize>>,
    pub callable_sites_per_pop_not_in_roh: Option<Array1<usize>>,

    // Per-population het tracking (summed across samples in pop)
    pub het_counts_per_pop: Array1<usize>,
    pub het_counts_per_pop_not_in_roh: Option<Array1<usize>>,

    // Per-site statistics (accumulated across all variant sites)
    pub variant_differences: Array1<usize>, // Pi differences per population
    pub variant_comparisons: Array1<usize>, // Pi comparisons per population

    // Dxy statistics (between populations)
    pub dxy_differences: Vec<usize>, // Between-pop differences
    pub dxy_comparisons: Vec<usize>, // Between-pop comparisons

    // Invariant comparisons
    pub invariant_comparisons: Array1<usize>,
    pub invariant_dxy_comparisons: Vec<usize>,

    pub fst_numerator: Vec<f64>,
    pub fst_denominator: Vec<f64>,

    // Track whether we have per-sample callable masks
    has_sample_masks: bool,
}

impl Window {
    pub fn new(
        region: Region,
        pop_map: PopulationMap,
        has_roh_data: bool,
        num_samples: usize,
        has_sample_masks: bool,
    ) -> Self {
        let num_pops = pop_map.num_populations();

        // Calculate number of population pairs for dxy
        let num_pop_pairs = if num_pops >= 2 {
            num_pops * (num_pops - 1) / 2
        } else {
            0
        };

        Self {
            region,
            pop_map,
            variant_positions: HashSet::new(),

            // Per-sample het tracking
            het_counts: Array1::zeros(num_samples),
            het_counts_not_in_roh: if has_roh_data {
                Some(Array1::zeros(num_samples))
            } else {
                None
            },

            // Per-sample callable sites (only for sample masks mode)
            callable_sites: if has_sample_masks {
                Some(Array1::zeros(num_samples))
            } else {
                None
            },
            callable_sites_not_in_roh: if has_sample_masks && has_roh_data {
                Some(Array1::zeros(num_samples))
            } else {
                None
            },

            // Per-population callable sites (for population counts mode)
            callable_sites_per_pop: if !has_sample_masks {
                Some(Array1::zeros(num_pops))
            } else {
                None
            },
            callable_sites_per_pop_not_in_roh: if !has_sample_masks && has_roh_data {
                Some(Array1::zeros(num_pops))
            } else {
                None
            },

            // Per-population het tracking
            het_counts_per_pop: Array1::zeros(num_pops),
            het_counts_per_pop_not_in_roh: if has_roh_data {
                Some(Array1::zeros(num_pops))
            } else {
                None
            },

            variant_differences: Array1::zeros(num_pops),
            variant_comparisons: Array1::zeros(num_pops),
            dxy_differences: vec![0; num_pop_pairs],
            dxy_comparisons: vec![0; num_pop_pairs],
            invariant_comparisons: Array1::zeros(num_pops),
            invariant_dxy_comparisons: vec![0; num_pop_pairs],
            fst_numerator: vec![0 as f64; num_pop_pairs],
            fst_denominator: vec![0 as f64; num_pop_pairs],
            has_sample_masks,
        }
    }

    /// Get index for population pair in the flat arrays
    fn get_pair_index(&self, pop1: usize, pop2: usize) -> usize {
        let num_pops = self.pop_map.num_populations();
        if pop1 < pop2 {
            pop1 * (num_pops - 1) - (pop1 * (pop1 + 1) / 2) + pop2 - 1
        } else {
            pop2 * (num_pops - 1) - (pop2 * (pop2 + 1) / 2) + pop1 - 1
        }
    }

    pub fn add_variant(&mut self, pos_idx: usize, counts: AlleleCounts) {
        // Track per-sample het counts (all samples)
        for (sample_idx, &is_het) in counts.hets.iter().enumerate() {
            if is_het {
                self.het_counts[sample_idx] += 1;
            }
        }

        // Track per-sample het counts (non-ROH only)
        if let Some(pos_non_roh) = &counts.non_roh {
            if let Some(het_counts_not_in_roh) = &mut self.het_counts_not_in_roh {
                for (sample_idx, &is_het) in pos_non_roh.hets.iter().enumerate() {
                    if is_het {
                        het_counts_not_in_roh[sample_idx] += 1;
                    }
                }
            }
        }

        // Track per-population het counts (summed)
        for pop_idx in 0..self.pop_map.num_populations() {
            // Count hets in this population
            let het_count_in_pop: usize = self
                .pop_map
                .get_population_by_index(pop_idx)
                .map(|pop| {
                    pop.samples()
                        .iter()
                        .filter_map(|sample_name| {
                            self.pop_map.lookup(sample_name).map(|(_, sample_idx)| {
                                if counts.hets.get(sample_idx).copied().unwrap_or(false) {
                                    1
                                } else {
                                    0
                                }
                            })
                        })
                        .sum()
                })
                .unwrap_or(0);
            self.het_counts_per_pop[pop_idx] += het_count_in_pop;

            // Non-ROH het counts per pop
            if let (Some(pos_non_roh), Some(het_counts_pop_not_in_roh)) =
                (&counts.non_roh, &mut self.het_counts_per_pop_not_in_roh)
            {
                let het_count_in_pop_non_roh: usize = self
                    .pop_map
                    .get_population_by_index(pop_idx)
                    .map(|pop| {
                        pop.samples()
                            .iter()
                            .filter_map(|sample_name| {
                                self.pop_map.lookup(sample_name).map(|(_, sample_idx)| {
                                    if pos_non_roh.hets.get(sample_idx).copied().unwrap_or(false) {
                                        1
                                    } else {
                                        0
                                    }
                                })
                            })
                            .sum()
                    })
                    .unwrap_or(0);
                het_counts_pop_not_in_roh[pop_idx] += het_count_in_pop_non_roh;
            }
        }

        // Track per-sample callable sites at variant positions
        if let Some(ref mut callable) = self.callable_sites {
            for (sample_idx, &is_callable) in counts.callable.iter().enumerate() {
                if is_callable {
                    callable[sample_idx] += 1;
                }
            }
        }

        // Track per-sample callable sites not in ROH at variant positions
        if let Some(pos_non_roh) = &counts.non_roh {
            if let Some(ref mut callable_not_in_roh) = self.callable_sites_not_in_roh {
                for (sample_idx, &is_callable) in pos_non_roh.callable.iter().enumerate() {
                    if is_callable {
                        callable_not_in_roh[sample_idx] += 1;
                    }
                }
            }
        }

        // Track per-population callable sites at variant positions (for population counts mode)
        if let Some(ref mut callable_per_pop) = self.callable_sites_per_pop {
            for pop_idx in 0..self.pop_map.num_populations() {
                let callable_count_in_pop: usize = self
                    .pop_map
                    .get_population_by_index(pop_idx)
                    .map(|pop| {
                        pop.samples()
                            .iter()
                            .filter_map(|sample_name| {
                                self.pop_map.lookup(sample_name).map(|(_, sample_idx)| {
                                    if counts.callable.get(sample_idx).copied().unwrap_or(false) {
                                        1
                                    } else {
                                        0
                                    }
                                })
                            })
                            .sum()
                    })
                    .unwrap_or(0);
                callable_per_pop[pop_idx] += callable_count_in_pop;
            }
        }

        // Track per-population callable sites not in ROH at variant positions
        if let Some(pos_non_roh) = &counts.non_roh {
            if let Some(ref mut callable_per_pop_not_in_roh) =
                self.callable_sites_per_pop_not_in_roh
            {
                for pop_idx in 0..self.pop_map.num_populations() {
                    let callable_count_in_pop: usize = self
                        .pop_map
                        .get_population_by_index(pop_idx)
                        .map(|pop| {
                            pop.samples()
                                .iter()
                                .filter_map(|sample_name| {
                                    self.pop_map.lookup(sample_name).map(|(_, sample_idx)| {
                                        if pos_non_roh
                                            .callable
                                            .get(sample_idx)
                                            .copied()
                                            .unwrap_or(false)
                                        {
                                            1
                                        } else {
                                            0
                                        }
                                    })
                                })
                                .sum()
                        })
                        .unwrap_or(0);
                    callable_per_pop_not_in_roh[pop_idx] += callable_count_in_pop;
                }
            }
        }

        let num_pops = counts.refs.len();
        let mut within_pi: Vec<f64> = Vec::with_capacity(num_pops);
        // Calculate per-site within-population statistics (Pi)
        for pop_idx in 0..num_pops {
            let refs = counts.refs[pop_idx];
            let alts = counts.alts[pop_idx];
            let total = refs + alts;

            // Per-site calculations
            let diffs = refs * alts;
            let comps = if total > 1 {
                (total * (total - 1)) / 2
            } else {
                0
            };

            self.variant_differences[pop_idx] += diffs;
            self.variant_comparisons[pop_idx] += comps;

            let pi = if comps > 0 {
                diffs as f64 / comps as f64
            } else {
                f64::NAN
            };
            within_pi.push(pi);
        }

        // Calculate per-site between-population statistics (Dxy)
        if num_pops >= 2 {
            for i in 0..num_pops {
                for j in (i + 1)..num_pops {
                    let refs1 = counts.refs[i];
                    let alts1 = counts.alts[i];
                    let refs2 = counts.refs[j];
                    let alts2 = counts.alts[j];

                    let total1 = refs1 + alts1;
                    let total2 = refs2 + alts2;

                    let diffs = refs1 * alts2 + alts1 * refs2;
                    let comps = total1 * total2;

                    let pair_idx = self.get_pair_index(i, j);
                    self.dxy_differences[pair_idx] += diffs;
                    self.dxy_comparisons[pair_idx] += comps;

                    let between = if comps > 0 {
                        diffs as f64 / comps as f64
                    } else {
                        f64::NAN
                    };

                    // Accumulate FST numerator/denominator per site
                    let pi_avg = (within_pi[i] + within_pi[j]) / 2.0;

                    if !between.is_nan() && !pi_avg.is_nan() && between > 0.0 {
                        let numerator = between - pi_avg;
                        let denominator = between;

                        self.fst_numerator[pair_idx] += numerator;
                        self.fst_denominator[pair_idx] += denominator;
                    }
                }
            }
        }

        self.variant_positions.insert(pos_idx);
    }

    pub fn add_invariant_comparisons(
        &mut self,
        comps: Array1<usize>,
        dxy_comps: Vec<usize>,
    ) {
        self.invariant_comparisons += &comps;

        for (idx, comp) in dxy_comps.iter().enumerate() {
            self.invariant_dxy_comparisons[idx] += comp;
        }
    }

    /// Add per-sample callable site counts from invariant positions (sample masks mode)
    pub fn add_callable_sites(
        &mut self,
        callable: Array1<usize>,
        callable_not_in_roh: Option<Array1<usize>>,
    ) {
        if let Some(ref mut sites) = self.callable_sites {
            *sites += &callable;
        }

        if let (Some(ref mut sites_not_in_roh), Some(callable_nr)) =
            (&mut self.callable_sites_not_in_roh, callable_not_in_roh)
        {
            *sites_not_in_roh += &callable_nr;
        }
    }

    /// Add per-population callable site counts from invariant positions (population counts mode)
    pub fn add_callable_sites_per_pop(
        &mut self,
        callable: Array1<usize>,
        callable_not_in_roh: Option<Array1<usize>>,
    ) {
        if let Some(ref mut sites) = self.callable_sites_per_pop {
            *sites += &callable;
        }

        if let (Some(ref mut sites_not_in_roh), Some(callable_nr)) =
            (&mut self.callable_sites_per_pop_not_in_roh, callable_not_in_roh)
        {
            *sites_not_in_roh += &callable_nr;
        }
    }
}

impl Window {
    fn compute_pi(&self) -> Vec<PiRecord> {
        let mut results = Vec::new();

        for (pop_idx, pop) in self.pop_map.populations().enumerate() {
            let differences = self.variant_differences[pop_idx];
            let comparisons =
                self.variant_comparisons[pop_idx] + self.invariant_comparisons[pop_idx];

            let pi = if comparisons > 0 {
                differences as f64 / comparisons as f64
            } else {
                f64::NAN
            };

            results.push(PiRecord::new(
                &self.region,
                pop.name.clone(),
                pi,
                comparisons,
                differences,
            ));
        }

        results
    }



    fn compute_dxy(&self) -> Option<Vec<DxyRecord>> {
        let num_pops = self.pop_map.num_populations();

        if num_pops < 2 {
            return None;
        }

        let mut results = Vec::new();

        for i in 0..num_pops {
            for j in (i + 1)..num_pops {
                let pop1 = self.pop_map.get_population_by_index(i).unwrap();
                let pop2 = self.pop_map.get_population_by_index(j).unwrap();

                let pair_idx = self.get_pair_index(i, j);
                let differences = self.dxy_differences[pair_idx];
                let comparisons =
                    self.dxy_comparisons[pair_idx] + self.invariant_dxy_comparisons[pair_idx];

                let dxy = if comparisons > 0 {
                    differences as f64 / comparisons as f64
                } else {
                    f64::NAN
                };

                results.push(DxyRecord::new(
                    &self.region,
                    pop1.name.clone(),
                    pop2.name.clone(),
                    dxy,
                    comparisons,
                    differences,
                ));
            }
        }

        Some(results)
    }

    fn compute_fst(&self) -> Option<Vec<FstRecord>> {
        let num_pops = self.pop_map.num_populations();

        if num_pops < 2 {
            return None;
        }

        let mut results = Vec::new();

        for i in 0..num_pops {
            for j in (i + 1)..num_pops {
                let pop1 = self.pop_map.get_population_by_index(i).unwrap();
                let pop2 = self.pop_map.get_population_by_index(j).unwrap();

                let pair_idx = self.get_pair_index(i, j);

                // FST = sum(numerators) / sum(denominators)
                let fst = if self.fst_denominator[pair_idx] > 0.0 {
                    self.fst_numerator[pair_idx] / self.fst_denominator[pair_idx]
                } else {
                    f64::NAN
                };

                results.push(FstRecord::new(
                    &self.region,
                    pop1.name.clone(),
                    pop2.name.clone(),
                    fst,
                ));
            }
        }

        Some(results)
    }
}

/// Heterozygosity record for output
#[derive(Debug, Serialize)]
pub struct HeterozygosityRecord {
    pub chrom: String,
    pub start: usize,
    pub end: usize,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub sample: Option<String>,
    pub population: String,
    pub het_total: usize,
    pub callable_total: usize,
    pub heterozygosity: f64,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub het_not_in_roh: Option<usize>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub callable_not_in_roh: Option<usize>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub heterozygosity_not_in_roh: Option<f64>,
}

impl HeterozygosityRecord {
    fn new_per_sample(
        region: &Region,
        sample: String,
        population: String,
        het_total: usize,
        callable_total: usize,
        het_not_in_roh: Option<usize>,
        callable_not_in_roh: Option<usize>,
    ) -> Self {
        let heterozygosity = if callable_total > 0 {
            het_total as f64 / callable_total as f64
        } else {
            f64::NAN
        };

        let heterozygosity_not_in_roh = match (het_not_in_roh, callable_not_in_roh) {
            (Some(het), Some(callable)) if callable > 0 => Some(het as f64 / callable as f64),
            (Some(_), Some(_)) => Some(f64::NAN),
            _ => None,
        };

        Self {
            chrom: region.name_as_str().to_string(),
            start: region.start_usize(),
            end: region.end_usize(),
            sample: Some(sample),
            population,
            het_total,
            callable_total,
            heterozygosity,
            het_not_in_roh,
            callable_not_in_roh,
            heterozygosity_not_in_roh,
        }
    }

    fn new_per_population(
        region: &Region,
        population: String,
        het_total: usize,
        callable_total: usize,
        het_not_in_roh: Option<usize>,
        callable_not_in_roh: Option<usize>,
    ) -> Self {
        let heterozygosity = if callable_total > 0 {
            het_total as f64 / callable_total as f64
        } else {
            f64::NAN
        };

        let heterozygosity_not_in_roh = match (het_not_in_roh, callable_not_in_roh) {
            (Some(het), Some(callable)) if callable > 0 => Some(het as f64 / callable as f64),
            (Some(_), Some(_)) => Some(f64::NAN),
            _ => None,
        };

        Self {
            chrom: region.name_as_str().to_string(),
            start: region.start_usize(),
            end: region.end_usize(),
            sample: None,
            population,
            het_total,
            callable_total,
            heterozygosity,
            het_not_in_roh,
            callable_not_in_roh,
            heterozygosity_not_in_roh,
        }
    }
}

#[derive(Debug, Serialize)]
pub struct PiRecord {
    pub chrom: String,
    pub start: usize,
    pub end: usize,
    population: String,
    pi: f64,
    comparisons: usize,
    differences: usize,
}
#[derive(Debug, Serialize)]
pub struct DxyRecord {
    pub chrom: String,
    pub start: usize,
    pub end: usize,
    population1: String,
    population2: String,
    dxy: f64,
    comparisons: usize,
    differences: usize,
}

#[derive(Debug, Serialize)]
pub struct FstRecord {
    pub chrom: String,
    pub start: usize,
    pub end: usize,
    population1: String,
    population2: String,
    fst: f64,
}

impl PiRecord {
    fn new(
        region: &Region,
        population: String,
        pi: f64,
        comparisons: usize,
        differences: usize,
    ) -> Self {
        Self {
            chrom: region.name_as_str().to_string(),
            start: region.start_usize(),
            end: region.end_usize(),
            population,
            pi,
            comparisons,
            differences,
        }
    }
}

impl DxyRecord {
    fn new(
        region: &Region,
        population1: String,
        population2: String,
        dxy: f64,
        comparisons: usize,
        differences: usize,
    ) -> Self {
        Self {
            chrom: region.name_as_str().to_string(),
            start: region.start_usize(),
            end: region.end_usize(),
            population1,
            population2,
            dxy,
            comparisons,
            differences,
        }
    }
}

impl FstRecord {
    fn new(region: &Region, population1: String, population2: String, fst: f64) -> Self {
        Self {
            chrom: region.name_as_str().to_string(),
            start: region.start_usize(),
            end: region.end_usize(),
            population1,
            population2,
            fst,
        }
    }
}

pub struct WindowStats {
    pub chrom: String,
    pub start: usize,
    pub end: usize,
    pub pi: Vec<PiRecord>,
    pub dxy: Option<Vec<DxyRecord>>,
    pub fst: Option<Vec<FstRecord>>,
    pub heterozygosity: Option<Vec<HeterozygosityRecord>>,
}

impl Window {
    pub fn compute_stats(&self, sample_names: &[String]) -> Result<WindowStats> {
        let pi = self.compute_pi();
        let dxy = self.compute_dxy();
        let fst = self.compute_fst();
        let heterozygosity = self.compute_heterozygosity(sample_names);
        Ok(WindowStats {
            chrom: self.region.name_as_str().to_string(),
            start: self.region.start_usize(),
            end: self.region.end_usize(),
            pi,
            dxy,
            fst,
            heterozygosity,
        })
    }

    /// Compute heterozygosity records
    fn compute_heterozygosity(&self, sample_names: &[String]) -> Option<Vec<HeterozygosityRecord>> {
        // Check if we have callable data (either per-sample or per-population)
        if self.callable_sites.is_none() && self.callable_sites_per_pop.is_none() {
            return None;
        }

        let mut results = Vec::new();

        if self.has_sample_masks {
            // Per-sample mode: output one record per sample
            if let Some(ref callable) = self.callable_sites {
                for (sample_idx, sample_name) in sample_names.iter().enumerate() {
                    // Find which population this sample belongs to
                    let pop_name = self
                        .pop_map
                        .lookup(sample_name)
                        .map(|(pop_idx, _)| {
                            self.pop_map
                                .get_population_by_index(pop_idx)
                                .map(|p| p.name.clone())
                                .unwrap_or_else(|| "unknown".to_string())
                        })
                        .unwrap_or_else(|| "unknown".to_string());

                    let het_total = self.het_counts[sample_idx];
                    let callable_total = callable[sample_idx];

                    let het_not_in_roh = self
                        .het_counts_not_in_roh
                        .as_ref()
                        .map(|h| h[sample_idx]);
                    let callable_not_in_roh = self
                        .callable_sites_not_in_roh
                        .as_ref()
                        .map(|c| c[sample_idx]);

                    results.push(HeterozygosityRecord::new_per_sample(
                        &self.region,
                        sample_name.clone(),
                        pop_name,
                        het_total,
                        callable_total,
                        het_not_in_roh,
                        callable_not_in_roh,
                    ));
                }
            }
        } else {
            // Per-population mode: output one record per population (summed approach)
            if let Some(ref callable) = self.callable_sites_per_pop {
                for (pop_idx, pop) in self.pop_map.populations().enumerate() {
                    let het_total = self.het_counts_per_pop[pop_idx];
                    let callable_total = callable[pop_idx];

                    let het_not_in_roh = self
                        .het_counts_per_pop_not_in_roh
                        .as_ref()
                        .map(|h| h[pop_idx]);
                    let callable_not_in_roh = self
                        .callable_sites_per_pop_not_in_roh
                        .as_ref()
                        .map(|c| c[pop_idx]);

                    results.push(HeterozygosityRecord::new_per_population(
                        &self.region,
                        pop.name.clone(),
                        het_total,
                        callable_total,
                        het_not_in_roh,
                        callable_not_in_roh,
                    ));
                }
            }
        }

        if results.is_empty() {
            None
        } else {
            Some(results)
        }
    }
}
