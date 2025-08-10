use super::regions::{CallableRegion, ChromRegion};

use bitvec::prelude::*;
use core::num;
use std::collections::HashMap;



#[derive(Debug, Clone)]
pub struct ChromosomeCounts {
    pub counts: BitVec,
    pub depth_sums: Vec<u32>,
    chrom: String,
    begin: u32,
    end: u32,
}

impl ChromosomeCounts {
    pub fn new(chrom: String, begin: u32, end: u32) -> Self {
        let len = (end - begin) as usize;
        Self {
            counts: bitvec![0; len],
            depth_sums: vec![0; len],
            chrom,
            begin,
            end,
        }
    }
}

pub struct GlobalCounts {
    chrom_counts: HashMap<String, AccumulatedCounts>,
}

pub struct GlobalPopCounts {
    chrom_counts: HashMap<String, AccumulatedPopCounts>,
    num_pops: u32

}



pub struct AccumulatedPopCounts {
    counts: Vec<Vec<u16>>,      // [position][population]
    depth_sums: Vec<Vec<u32>>,  // [position][population]
    chrom: String,
    begin: u32,
    end: u32,
    num_pops: usize,
}

impl AccumulatedPopCounts {
    pub fn new(chrom: String, begin: u32, end: u32, num_pops: usize) -> Self {
        let len = (end - begin) as usize;
        Self {
            counts: vec![vec![0; num_pops]; len],
            depth_sums: vec![vec![0; num_pops]; len],
            chrom,
            begin,
            end,
            num_pops,
        }
    }

    pub fn add_chromosome_counts(&mut self, counts: ChromosomeCounts, pop_idx: usize) {
        assert_eq!(self.chrom, counts.chrom);
        assert_eq!(self.begin, counts.begin);
        assert_eq!(self.end, counts.end);

        for (i, is_set) in counts.counts.iter().by_vals().enumerate() {
            if is_set {
                self.counts[i][pop_idx] = self.counts[i][pop_idx].saturating_add(1);
            }
        }
        for (acc_depths, &new_depth) in self.depth_sums.iter_mut().zip(counts.depth_sums.iter()) {
            acc_depths[pop_idx] = acc_depths[pop_idx].saturating_add(new_depth);
        }
    }

    
    
}

impl AccumulatedPopCounts {
    /// For each population, output a region for every base that passes the thresholds.
    pub fn to_callable_regions(
        &self,
        min_proportion: f64,
        mean_thresholds: (f64, f64),
        total_samples: usize,
    ) -> Vec<Vec<CallableRegion>> {
        let mut regions_per_pop = vec![Vec::new(); self.num_pops];
        let min_samples = (total_samples as f64 * min_proportion).ceil() as u32;

        for (i, (counts, depths)) in self.counts.iter().zip(self.depth_sums.iter()).enumerate() {
            let total_count: u32 = counts.iter().map(|&c| c as u32).sum();
            let total_depth: u32 = depths.iter().sum();
            let mean_depth = if total_count > 0 {
                total_depth as f64 / total_count as f64
            } else {
                0.0
            };

            if total_count >= min_samples
                && mean_depth >= mean_thresholds.0
                && mean_depth <= mean_thresholds.1
            {
                let pos = self.begin + i as u32;
                for (pop_idx, &count) in counts.iter().enumerate() {
                    regions_per_pop[pop_idx].push(CallableRegion {
                        begin: pos,
                        end: pos + 1,
                        count: count as u32,
                    });
                }
            }
        }

        regions_per_pop
    }
}
pub struct AccumulatedCounts {
    counts: Vec<u16>,
    depth_sums: Vec<u32>,
    chrom: String,
    begin: u32,
    end: u32,
    
}

impl AccumulatedCounts {
    fn new(chrom: String, begin: u32, end: u32) -> Self {
        let len = (end - begin) as usize;
        Self {
            counts: vec![0; len],
            depth_sums: vec![0; len],
            chrom,
            begin,
            end,
        }
    }

    fn add_chromosome_counts(&mut self, counts: ChromosomeCounts) {
        assert_eq!(self.chrom, counts.chrom);
        assert_eq!(self.begin, counts.begin);
        assert_eq!(self.end, counts.end);

        for (i, is_set) in counts.counts.iter().by_vals().enumerate() {
            if is_set {
                self.counts[i] = self.counts[i].saturating_add(1);
            }
        }

        for (acc_depth, &new_depth) in self.depth_sums.iter_mut().zip(counts.depth_sums.iter()) {
            *acc_depth = acc_depth.saturating_add(new_depth);
        }
    }

    fn to_callable_regions(
        &self,
        min_samples: u32,
        mean_thresolds: (f64, f64),
    ) -> Vec<CallableRegion> {
        let mut regions = Vec::new();
        let mut current_start = None;
        let mut current_stats = None;

        let mut add_region = |start: u32, end: u32, count: u16, depth_sum: u32| {
            let mean_depth = if count > 0 {
                depth_sum as f64 / count as f64
            } else {
                0.0
            };

            if count >= min_samples as u16
                && mean_depth >= mean_thresolds.0
                && mean_depth <= mean_thresolds.1
            {
                regions.push(CallableRegion {
                    begin: start,
                    end,
                    count: count as u32,
                });
            }
        };

        for (i, (&count, &depth_sum)) in self.counts.iter().zip(self.depth_sums.iter()).enumerate()
        {
            let pos = self.begin + i as u32;

            match (current_start, current_stats) {
                (Some(start), Some((prev_count, prev_depth)))
                    if count == prev_count && depth_sum == prev_depth =>
                {
                    continue;
                }
                (Some(start), Some((prev_count, prev_depth))) => {
                    add_region(start, pos, prev_count, prev_depth);
                    if count > 0 {
                        current_start = Some(pos);
                        current_stats = Some((count, depth_sum));
                    } else {
                        current_start = None;
                        current_stats = None;
                    }
                }
                (None, None) if count > 0 => {
                    current_start = Some(pos);
                    current_stats = Some((count, depth_sum));
                }
                _ => {}
            }
        }

        if let (Some(start), Some((count, depth_sum))) = (current_start, current_stats) {
            add_region(start, self.end, count, depth_sum);
        }

        regions
    }
}

impl GlobalCounts {
    pub fn new(regions: &[ChromRegion]) -> Self {
        let mut chrom_counts = HashMap::new();

        for region in regions {
            chrom_counts.insert(
                region.chr.clone(),
                AccumulatedCounts::new(region.chr.clone(), region.begin, region.end),
            );
        }

        Self { chrom_counts }
    }

    pub fn merge(&mut self, chromosome_counts: ChromosomeCounts) {
        if let Some(existing) = self.chrom_counts.get_mut(&chromosome_counts.chrom) {
            existing.add_chromosome_counts(chromosome_counts);
        } else {
            let mut accumulated = AccumulatedCounts::new(
                chromosome_counts.chrom.clone(),
                chromosome_counts.begin,
                chromosome_counts.end,
            );
            accumulated.add_chromosome_counts(chromosome_counts);
            self.chrom_counts
                .insert(accumulated.chrom.clone(), accumulated);
        }
    }

    pub fn finalize(
        self,
        min_proportion: f64,
        mean_thresholds: (f64, f64),
        total_samples: usize,
    ) -> Vec<(String, u32, Vec<CallableRegion>)> {
        let min_samples = (total_samples as f64 * min_proportion).ceil() as u32;

        let mut results: Vec<_> = self
            .chrom_counts
            .into_iter()
            .map(|(chrom, counts)| {
                let begin = counts.begin;
                let regions = counts.to_callable_regions(min_samples, mean_thresholds);
                (chrom, begin, regions)
            })
            .collect();

        results.sort_by(|a, b| a.0.cmp(&b.0));
        results
    }
}


impl GlobalPopCounts {
    pub fn new(regions: &[ChromRegion], num_pops:u32) -> Self {
        let mut chrom_counts = HashMap::new();

        for region in regions {
            chrom_counts.insert(
                region.chr.clone(),
                AccumulatedPopCounts::new(region.chr.clone(), region.begin, region.end, num_pops as usize),
            );
        }

        Self { chrom_counts, num_pops }
    }

    pub fn merge(&mut self, chromosome_counts: ChromosomeCounts, pop_idx:usize) {
        if let Some(existing) = self.chrom_counts.get_mut(&chromosome_counts.chrom) {
            existing.add_chromosome_counts(chromosome_counts, pop_idx);
        } else {
            let mut accumulated = AccumulatedPopCounts::new(
                chromosome_counts.chrom.clone(),
                chromosome_counts.begin,
                chromosome_counts.end,
                self.num_pops as usize,
            );
            accumulated.add_chromosome_counts(chromosome_counts, pop_idx);
            self.chrom_counts
                .insert(accumulated.chrom.clone(), accumulated);
        }
    }

    pub fn finalize(
    &self,
    min_proportion: f64,
    mean_thresholds: (f64, f64),
    total_samples: usize,
    pop_idx: usize,
) -> Vec<(String, u32, Vec<CallableRegion>)> {
    let mut results: Vec<_> = self
        .chrom_counts
        .iter()
        .map(|(chrom, counts)| {
            let begin = counts.begin;
            let regions_per_pop = counts.to_callable_regions(min_proportion, mean_thresholds, total_samples);
            // Only take the regions for the requested population
            (chrom.clone(), begin, regions_per_pop[pop_idx].clone())
        })
        .collect();

    results.sort_by(|a, b| a.0.cmp(&b.0));
    results
}
}
#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_chromosome_counts_new() {
        let counts = ChromosomeCounts::new("chr1".to_string(), 1000, 2000);
        assert_eq!(counts.counts.len(), 1000);
        assert_eq!(counts.depth_sums.len(), 1000);
        assert!(counts.counts.not_any()); // All bits should be 0
        assert!(counts.depth_sums.iter().all(|&x| x == 0));
    }

    #[test]
    fn test_accumulated_counts_new() {
        let acc = AccumulatedCounts::new("chr1".to_string(), 1000, 2000);
        assert_eq!(acc.counts.len(), 1000);
        assert_eq!(acc.depth_sums.len(), 1000);
        assert!(acc.counts.iter().all(|&x| x == 0));
        assert!(acc.depth_sums.iter().all(|&x| x == 0));
    }

    #[test]
    fn test_add_chromosome_counts() {
        let mut acc = AccumulatedCounts::new("chr1".to_string(), 1000, 2000);
        let mut chrom_counts = ChromosomeCounts::new("chr1".to_string(), 1000, 2000);

        // Set some bits and depths
        chrom_counts.counts.set(0, true);
        chrom_counts.counts.set(100, true);
        chrom_counts.depth_sums[0] = 5;
        chrom_counts.depth_sums[100] = 10;

        acc.add_chromosome_counts(chrom_counts);

        assert_eq!(acc.counts[0], 1);
        assert_eq!(acc.counts[100], 1);
        assert_eq!(acc.depth_sums[0], 5);
        assert_eq!(acc.depth_sums[100], 10);

        // Test counts don't exceed u16::MAX
        let mut max_counts = ChromosomeCounts::new("chr1".to_string(), 1000, 2000);
        max_counts.counts.set(0, true);
        max_counts.depth_sums[0] = 1;

        for _ in 0..=u16::MAX as usize {
            acc.add_chromosome_counts(max_counts.clone());
        }

        assert_eq!(acc.counts[0], u16::MAX);
    }

    #[test]
    fn test_global_counts_merge() {
        let regions = vec![ChromRegion {
            chr: "chr1".to_string(),
            begin: 1000,
            end: 2000,
            min_filter: 0.0,
            max_filter: 0.0,
        }];
        let mut global = GlobalCounts::new(&regions);

        let mut counts1 = ChromosomeCounts::new("chr1".to_string(), 1000, 2000);
        counts1.counts.set(50, true); // Position 1050
        counts1.depth_sums[50] = 5;

        let mut counts2 = ChromosomeCounts::new("chr1".to_string(), 1000, 2000);
        counts2.counts.set(50, true); // Same position
        counts2.depth_sums[50] = 3;

        global.merge(counts1);
        global.merge(counts2);

        let results = global.finalize(0.1, (1.0, f64::MAX), 2); // min_proportion = 0.1 means min_samples = 1
        assert_eq!(results.len(), 1); // One chromosome

        let (_chrom, _begin, regions) = &results[0];
        assert!(!regions.is_empty());

        // Find region containing our test position
        let region = regions.iter().find(|r| r.begin <= 1050 && r.end > 1050);
        assert!(region.is_some());
        let region = region.unwrap();
        assert_eq!(region.count, 2); // Both samples contributed
    }

    #[test]
    fn test_callable_regions_filtering() {
        let mut acc = AccumulatedCounts::new("chr1".to_string(), 1000, 2000);

        // Create regions with different counts and depths
        acc.counts[0] = 5; // Position 1000
        acc.depth_sums[0] = 25; // Mean depth 5

        acc.counts[100] = 5; // Position 1100
        acc.depth_sums[100] = 10; // Mean depth 2

        acc.counts[200] = 2; // Position 1200
        acc.depth_sums[200] = 10; // Mean depth 5

        // Test with different thresholds
        let regions1 = acc.to_callable_regions(3, (4.0, f64::MAX)); // Should include pos 1000 only
        let regions2 = acc.to_callable_regions(3, (1.0, f64::MAX)); // Should include pos 1000 and 1100

        assert_eq!(regions1.len(), 1);
        assert_eq!(regions2.len(), 2);
    }
}
