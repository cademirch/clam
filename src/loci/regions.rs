use std::collections::{HashMap, HashSet};

use anyhow::{bail, Result};
use log::warn;

use super::thresholds::Thresholds;

#[derive(Clone)]
pub struct ChromRegion {
    pub chr: String,
    pub begin: u32,
    pub end: u32,
    pub min_filter: f64,
    pub max_filter: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub struct CallableRegion {
    pub count: u32,
    pub begin: u32,
    pub end: u32,
}

pub fn intersect_intervals(
    pop_regions: &[CallableRegion],
    global_regions: &[CallableRegion],
) -> Vec<CallableRegion> {
    let mut result = Vec::new();
    let mut i = 0;
    let mut j = 0;

    while i < pop_regions.len() && j < global_regions.len() {
        let a = &pop_regions[i];
        let b = &global_regions[j];

        let start = a.begin.max(b.begin);
        let end = a.end.min(b.end);

        if start < end {
            // There is an intersection
            result.push(CallableRegion {
                begin: start,
                end,
                count: a.count, // or min(a.count, b.count) if you want
            });
        }

        // Move to next interval
        if a.end < b.end {
            i += 1;
        } else {
            j += 1;
        }
    }
    result
}

pub fn prepare_chrom_regions(
    chrom_regions: Vec<(&str, u32, u32)>,
    thresholds: Thresholds,
    exclude_chrs: Option<&HashSet<String>>,
    include_chrs: Option<&HashSet<String>>,
) -> Result<Vec<ChromRegion>> {
    // Validate that all included chromosomes exist in chrom_regions
    if let Some(include_set) = include_chrs {
        for chrom in include_set.iter() {
            if !chrom_regions
                .iter()
                .any(|(region_chrom, _, _)| *region_chrom == *chrom)
            {
                bail!(
                    "Included chromosome '{}' not found in input file's chrom regions.",
                    chrom
                );
            }
        }
    }

    let regions: Vec<ChromRegion> = match thresholds {
        Thresholds::Fixed((min_depth, max_depth)) => chrom_regions
            .into_iter()
            .map(|(chr, start, end)| ChromRegion {
                chr: chr.to_string(),
                begin: start,
                end,
                min_filter: min_depth,
                max_filter: max_depth,
            })
            .collect(),

        Thresholds::PerChromosome(filter_map) => {
            // First, validate that all filters have corresponding chromosomes
            for (chrom, _) in filter_map.iter() {
                if !chrom_regions
                    .iter()
                    .any(|(region_chrom, _, _)| region_chrom == chrom)
                {
                    bail!(
                        "Filter provided for chromosome '{}' not found in D4 chrom regions.",
                        chrom
                    );
                }
            }

            // Then create ChromRegions with appropriate filters
            chrom_regions
                .into_iter()
                .map(|(chrom, start, end)| {
                    let (min_filter, max_filter) = filter_map
                        .get(chrom)
                        .copied()
                        .unwrap_or_else(|| {
                            warn!(
                                "No filter found for chromosome '{}', using default filters (0.0, f64::INFINITY).",
                                chrom
                            );
                            (0.0, f64::INFINITY)
                        });

                    ChromRegion {
                        chr: chrom.to_string(),
                        begin: start,
                        end,
                        min_filter,
                        max_filter,
                    }
                })
                .collect()
        }
    };

    // Filter based on both include and exclude lists
    let filtered_regions: Vec<ChromRegion> = regions
        .into_iter()
        .filter(|region| {
            let is_not_excluded =
                exclude_chrs.map_or(true, |exclude_set| !exclude_set.contains(&region.chr));
            let is_included =
                include_chrs.map_or(true, |include_set| include_set.contains(&region.chr));
            is_not_excluded && is_included
        })
        .collect();

    // Check if any regions remain after filtering
    if filtered_regions.is_empty() {
        bail!("No chromosome regions remaining after applying include/exclude filters");
    }

    Ok(filtered_regions)
}
