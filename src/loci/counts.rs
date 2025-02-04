use super::CallableRegion;
use super::ChromRegion;
use std::collections::HashMap;

/// Stores counts for each position in a chromosome
pub struct ChromosomeCounts {
    counts: Vec<u32>,
    chrom: String,
    begin: u32,
    end: u32,
}

impl ChromosomeCounts {
    pub fn new(chrom: String, begin: u32, end: u32) -> Self {
        let len = (end - begin) as usize;
        Self {
            counts: vec![0; len],
            chrom,
            begin,
            end,
        }
    }

    pub fn add_regions(&mut self, regions: &[CallableRegion]) {
        for region in regions {
            if region.begin >= self.end || region.end <= self.begin {
                continue;
            }
            
            let start_idx = (region.begin.saturating_sub(self.begin)) as usize;
            let end_idx = ((region.end - self.begin).min(self.counts.len() as u32)) as usize;
            
            for count in &mut self.counts[start_idx..end_idx] {
                *count += 1;
            }
        }
    }

    pub fn to_callable_regions(
        &self,
        min_proportion: f64,
        total_samples: usize,
        output_counts: bool,
    ) -> Vec<CallableRegion> {
        let mut regions = Vec::new();
        let min_samples = if output_counts { 
            1 
        } else {
            (total_samples as f64 * min_proportion).ceil() as u32
        };

        let mut current_start = None;
        let mut current_count = None;

        // Helper to add a region to our results
        let mut add_region = |start: u32, end: u32, count: u32| {
            if count >= min_samples || output_counts {
                regions.push(CallableRegion {
                    begin: start,
                    end,
                    count,
                });
            }
        };

        for (i, &count) in self.counts.iter().enumerate() {
            let pos = self.begin + i as u32;

            match (current_start, current_count) {
                (Some(start), Some(prev_count)) if count == prev_count => {
                    // Continue current region
                    continue;
                }
                (Some(start), Some(prev_count)) => {
                    // End current region and maybe start new one
                    add_region(start, pos, prev_count);
                    if count > 0 {
                        current_start = Some(pos);
                        current_count = Some(count);
                    } else {
                        current_start = None;
                        current_count = None;
                    }
                }
                (None, None) if count > 0 => {
                    // Start new region
                    current_start = Some(pos);
                    current_count = Some(count);
                }
                _ => {}
            }
        }

        // Handle last region if exists
        if let (Some(start), Some(count)) = (current_start, current_count) {
            add_region(start, self.end, count);
        }

        regions
    }
}

/// Global state for tracking counts across all chromosomes
pub struct GlobalCounts {
    chrom_counts: HashMap<String, ChromosomeCounts>,
    total_samples: usize,
}

impl GlobalCounts {
    pub fn new(regions: &[ChromRegion]) -> Self {
        let mut chrom_counts = HashMap::new();
        
        for region in regions {
            chrom_counts.insert(
                region.chr.clone(),
                ChromosomeCounts::new(region.chr.clone(), region.begin, region.end)
            );
        }
        
        Self {
            chrom_counts,
            total_samples: 0,
        }
    }

    pub fn add_sample_regions(&mut self, sample_regions: Vec<(String, Vec<CallableRegion>)>) {
        self.total_samples += 1;
        
        for (chrom, regions) in sample_regions {
            if let Some(counts) = self.chrom_counts.get_mut(&chrom) {
                counts.add_regions(&regions);
            }
        }
    }

    pub fn finalize(
        self,
        min_proportion: f64,
        output_counts: bool,
    ) -> Vec<(String, u32, Vec<CallableRegion>)> {
        let mut results: Vec<_> = self.chrom_counts
            .into_iter()
            .map(|(chrom, counts)| {
                let regions = counts.to_callable_regions(
                    min_proportion,
                    self.total_samples,
                    output_counts,
                );
                (chrom, 0, regions)
            })
            .collect();
            
        results.sort_by(|a, b| a.0.cmp(&b.0));
        results
    }
}