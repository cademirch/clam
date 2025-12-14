use crate::core::population::PopulationMap;
use crate::stat::vcf::variants::AlleleCounts;
use bstr::ByteSlice;
use color_eyre::Result;
use color_eyre::eyre::Ok;
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
    pub total_counts: AlleleCounts,
    pub variant_positions: HashSet<usize>,
    pub het_counts: Array1<usize>,
    pub het_counts_non_roh: Option<Array1<usize>>,
}

impl Window {
    pub fn new(region: Region, pop_map: PopulationMap, has_roh_data: bool) -> Self {
        let num_pops = pop_map.num_populations();
        let num_samples = pop_map.all_samples().len();

        // Initialize empty AlleleCounts
        let total_counts = AlleleCounts {
            refs: Array1::zeros(num_pops),
            alts: Array1::zeros(num_pops),
            missing: Array1::zeros(num_pops),
            hets: Array1::from_elem(num_samples, false),
            non_roh: if has_roh_data {
                Some(crate::stat::vcf::variants::NonRohCounts {
                    refs: Array1::zeros(num_pops),
                    alts: Array1::zeros(num_pops),
                    hets: Array1::from_elem(num_samples, false),
                })
            } else {
                None
            },
        };

        Self {
            region,
            pop_map,
            total_counts,
            variant_positions: HashSet::new(),
            het_counts: Array1::zeros(num_samples),
            het_counts_non_roh: if has_roh_data {
                Some(Array1::zeros(num_samples))
            } else {
                None
            },
        }
    }

    pub fn add_variant(&mut self, pos_idx: usize, counts: AlleleCounts) {
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

        self.total_counts += counts;

        self.variant_positions.insert(pos_idx);
    }

    pub fn add_invariants(&mut self, refs_all: Array1<usize>, refs_non_roh: Option<Array1<usize>>) {
        self.total_counts.refs += &refs_all;

        if let (Some(non_roh), Some(refs_non_roh)) = (&mut self.total_counts.non_roh, refs_non_roh)
        {
            non_roh.refs += &refs_non_roh;
        }
    }
}

#[derive(Debug, Serialize)]
pub struct PiRecord {
    #[serde(flatten)]
    region: SerializableRegion,
    population: String,
    pi: f64,
    comparisons: usize,
    differences: usize,
}
#[derive(Debug, Serialize)]
pub struct DxyRecord {
    #[serde(flatten)]
    region: SerializableRegion,
    population1: String,
    population2: String,
    dxy: f64,
    comparisons: usize,
    differences: usize,
}

#[derive(Debug, Serialize)]
pub struct FstRecord {
    #[serde(flatten)]
    region: SerializableRegion,
    population1: String,
    population2: String,
    fst: f64,
}
impl Window {
    fn compute_pi(&self) -> Vec<PiRecord> {
        let mut results = Vec::new();

        for (pop_idx, pop) in self.pop_map.populations().enumerate() {
            let refs = self.total_counts.refs[pop_idx];
            let alts = self.total_counts.alts[pop_idx];
            let total = refs + alts;

            let (pi, comparisons, differences) = if total > 1 {
                let differences = refs * alts;
                let comparisons = total * (total - 1) / 2;
                let pi = differences as f64 / comparisons as f64;
                (pi, comparisons, differences)
            } else {
                (f64::NAN, 0, 0)
            };

            results.push(PiRecord {
                population: pop.name.clone(),
                region: SerializableRegion::from(&self.region),
                pi,
                comparisons,
                differences,
            });
        }

        results
    }
    fn compute_pi_non_roh(&self) -> Option<Vec<PiRecord>> {
        let results = if let Some(non_roh) = self.total_counts.non_roh.as_ref() {
            let mut results = Vec::new();

            for (pop_idx, pop) in self.pop_map.populations().enumerate() {
                let refs = non_roh.refs[pop_idx];
                let alts = non_roh.alts[pop_idx];
                let total = refs + alts;

                let (pi, comparisons, differences) = if total > 1 {
                    let differences = refs * alts;
                    let comparisons = total * (total - 1) / 2;
                    let pi = differences as f64 / comparisons as f64;
                    (pi, comparisons, differences)
                } else {
                    (f64::NAN, 0, 0)
                };

                results.push(PiRecord {
                    population: pop.name.clone(),
                    region: SerializableRegion::from(&self.region),
                    pi,
                    comparisons,
                    differences,
                });
            }

            Some(results)
        } else {
            None
        };

        results
    }

    fn compute_dxy(&self) -> Option<Vec<DxyRecord>> {
        let mut results = Vec::new();
        let num_pops = self.pop_map.num_populations();

        if num_pops < 2 {
            return None;
        }

        for i in 0..num_pops {
            for j in (i + 1)..num_pops {
                let pop1 = self.pop_map.get_population_by_index(i).unwrap();
                let pop2 = self.pop_map.get_population_by_index(j).unwrap();

                let refs1 = self.total_counts.refs[i];
                let alts1 = self.total_counts.alts[i];
                let refs2 = self.total_counts.refs[j];
                let alts2 = self.total_counts.alts[j];

                let total1 = refs1 + alts1;
                let total2 = refs2 + alts2;

                let (dxy, comparisons, differences) = if total1 > 0 && total2 > 0 {
                    let comparisons = total1 * total2;
                    let differences = refs1 * alts2 + alts1 * refs2;
                    let dxy = differences as f64 / comparisons as f64;
                    (dxy, comparisons, differences)
                } else {
                    (f64::NAN, 0, 0)
                };

                results.push(DxyRecord {
                    region: SerializableRegion::from(&self.region),
                    population1: pop1.name.clone(),
                    population2: pop2.name.clone(),
                    dxy,
                    comparisons,
                    differences,
                });
            }
        }

        Some(results)
    }

    fn compute_fst(&self) -> Option<Vec<FstRecord>> {
        let mut results = Vec::new();
        let num_pops = self.pop_map.num_populations();

        if num_pops < 2 {
            return None;
        }

        let mut pi_values = Vec::new();
        for pop_idx in 0..num_pops {
            let refs = self.total_counts.refs[pop_idx];
            let alts = self.total_counts.alts[pop_idx];
            let total = refs + alts;

            let pi = if total > 1 {
                let differences = refs * alts;
                let comparisons = total * (total - 1) / 2;
                differences as f64 / comparisons as f64
            } else {
                f64::NAN
            };
            pi_values.push(pi);
        }

        for i in 0..num_pops {
            for j in (i + 1)..num_pops {
                let pop1 = self.pop_map.get_population_by_index(i).unwrap();
                let pop2 = self.pop_map.get_population_by_index(j).unwrap();

                let refs1 = self.total_counts.refs[i];
                let alts1 = self.total_counts.alts[i];
                let refs2 = self.total_counts.refs[j];
                let alts2 = self.total_counts.alts[j];

                let total1 = refs1 + alts1;
                let total2 = refs2 + alts2;

                let fst = if total1 > 0 && total2 > 0 {
                    let comparisons = total1 * total2;
                    let differences = refs1 * alts2 + alts1 * refs2;
                    let dxy = differences as f64 / comparisons as f64;

                    let pi_avg = (pi_values[i] + pi_values[j]) / 2.0;

                    // FST = (dxy - pi_avg) / dxy
                    if dxy > 0.0 {
                        (dxy - pi_avg) / dxy
                    } else {
                        f64::NAN
                    }
                } else {
                    f64::NAN
                };

                results.push(FstRecord {
                    region: SerializableRegion::from(&self.region),
                    population1: pop1.name.clone(),
                    population2: pop2.name.clone(),
                    fst,
                });
            }
        }

        Some(results)
    }
}

pub struct WindowStats {
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
        Ok(WindowStats { pi, pi_non_roh, dxy, fst })
    }
}

