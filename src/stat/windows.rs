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

        let num_windows = ((end - start) + window_size - 1) / window_size;

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
        let end = (start + self.window_size).min(self.region_end);
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
    pub het_counts: Array1<usize>,
    pub het_counts_non_roh: Option<Array1<usize>>,

    // Per-site statistics (accumulated across all variant sites)
    pub variant_differences: Array1<usize>, // Pi differences per population
    pub variant_comparisons: Array1<usize>, // Pi comparisons per population
    pub variant_differences_non_roh: Option<Array1<usize>>,
    pub variant_comparisons_non_roh: Option<Array1<usize>>,

    // Dxy statistics (between populations)
    pub dxy_differences: Vec<usize>, // Between-pop differences
    pub dxy_comparisons: Vec<usize>, // Between-pop comparisons

    // Invariant comparisons
    pub invariant_comparisons: Array1<usize>,
    pub invariant_comparisons_non_roh: Option<Array1<usize>>,
    pub invariant_dxy_comparisons: Vec<usize>,

    pub fst_numerator: Vec<f64>,
    pub fst_denominator: Vec<f64>,
}

impl Window {
    pub fn new(region: Region, pop_map: PopulationMap, has_roh_data: bool) -> Self {
        let num_pops = pop_map.num_populations();
        let num_samples = pop_map.all_samples().len();

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
            het_counts: Array1::zeros(num_samples),
            het_counts_non_roh: if has_roh_data {
                Some(Array1::zeros(num_samples))
            } else {
                None
            },
            variant_differences: Array1::zeros(num_pops),
            variant_comparisons: Array1::zeros(num_pops),
            variant_differences_non_roh: if has_roh_data {
                Some(Array1::zeros(num_pops))
            } else {
                None
            },
            variant_comparisons_non_roh: if has_roh_data {
                Some(Array1::zeros(num_pops))
            } else {
                None
            },
            dxy_differences: vec![0; num_pop_pairs],
            dxy_comparisons: vec![0; num_pop_pairs],
            invariant_comparisons: Array1::zeros(num_pops),
            invariant_comparisons_non_roh: if has_roh_data {
                Some(Array1::zeros(num_pops))
            } else {
                None
            },
            invariant_dxy_comparisons: vec![0; num_pop_pairs],
            fst_numerator: vec![0 as f64; num_pop_pairs],
            fst_denominator: vec![0 as f64; num_pop_pairs],
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
        // Track het counts
        for (sample_idx, &is_het) in counts.hets.iter().enumerate() {
            if is_het {
                self.het_counts[sample_idx] += 1;
            }
        }

        if let Some(pos_non_roh) = &counts.non_roh {
            if let Some(het_counts_non_roh) = &mut self.het_counts_non_roh {
                for (sample_idx, &is_het) in pos_non_roh.hets.iter().enumerate() {
                    if is_het {
                        het_counts_non_roh[sample_idx] += 1;
                    }
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

        // Non-ROH calculations
        if let Some(non_roh) = &counts.non_roh {
            if let (Some(diffs_non_roh), Some(comps_non_roh)) = (
                &mut self.variant_differences_non_roh,
                &mut self.variant_comparisons_non_roh,
            ) {
                for pop_idx in 0..non_roh.refs.len() {
                    let refs = non_roh.refs[pop_idx];
                    let alts = non_roh.alts[pop_idx];
                    let total = refs + alts;

                    let diffs = refs * alts;
                    let comps = if total > 1 {
                        (total * (total - 1)) / 2
                    } else {
                        0
                    };

                    diffs_non_roh[pop_idx] += diffs;
                    comps_non_roh[pop_idx] += comps;
                }
            }
        }

        self.variant_positions.insert(pos_idx);
    }

    pub fn add_invariant_comparisons(
        &mut self,
        comps: Array1<usize>,
        comps_non_roh: Option<Array1<usize>>,
        dxy_comps: Vec<usize>,
    ) {
        self.invariant_comparisons += &comps;

        for (idx, comp) in dxy_comps.iter().enumerate() {
            self.invariant_dxy_comparisons[idx] += comp;
        }

        if let (Some(inv_comps_non_roh), Some(comps_non_roh)) =
            (&mut self.invariant_comparisons_non_roh, comps_non_roh)
        {
            *inv_comps_non_roh += &comps_non_roh;
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

    fn compute_pi_non_roh(&self) -> Option<Vec<PiRecord>> {
        if let (Some(diffs), Some(var_comps), Some(inv_comps)) = (
            &self.variant_differences_non_roh,
            &self.variant_comparisons_non_roh,
            &self.invariant_comparisons_non_roh,
        ) {
            let mut results = Vec::new();

            for (pop_idx, pop) in self.pop_map.populations().enumerate() {
                let differences = diffs[pop_idx];
                let comparisons = var_comps[pop_idx] + inv_comps[pop_idx];

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

            Some(results)
        } else {
            None
        }
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
    pub pi_non_roh: Option<Vec<PiRecord>>,
    pub dxy: Option<Vec<DxyRecord>>,
    pub fst: Option<Vec<FstRecord>>,
}

impl Window {
    pub fn compute_stats(&self) -> Result<WindowStats> {
        let pi = self.compute_pi();
        let dxy = self.compute_dxy();
        let fst = self.compute_fst();
        let pi_non_roh = self.compute_pi_non_roh();
        Ok(WindowStats {
            chrom: self.region.name_as_str().to_string(),
            start: self.region.start_usize(),
            end: self.region.end_usize(),
            pi,
            pi_non_roh,
            dxy,
            fst,
        })
    }
}
