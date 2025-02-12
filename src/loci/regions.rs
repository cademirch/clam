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

pub fn prepare_chrom_regions(
    chrom_regions: Vec<(&str, u32, u32)>,
    thresholds: Thresholds,
    exclude_chrs: Option<&HashSet<String>>,
) -> Result<Vec<ChromRegion>> {
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

    // Filter out excluded chromosomes
    Ok(regions
        .into_iter()
        .filter(|region| {
            exclude_chrs.map_or(true, |exclude_set| !exclude_set.contains(&region.chr))
        })
        .collect())
}
pub fn add_filters_to_chroms(
    chrom_regions: Vec<(&str, u32, u32)>,
    filter_map: HashMap<String, (f64, f64)>,
) -> Result<Vec<(String, u32, u32, f64, f64)>> {
    let mut result = vec![];
    for (chrom, start, end) in &chrom_regions {
        if let Some(&(min_filter, max_filter)) = filter_map.get(*chrom) {
            result.push((chrom.to_string(), *start, *end, min_filter, max_filter));
        } else {
            // Issue a warning if no filter is found and apply default filters
            warn!(
                "No filter found for chromosome '{}', using default filters (0.0, f64::INFINITY).",
                chrom
            );
            result.push((chrom.to_string(), *start, *end, 0.0, f64::INFINITY));
        }
    }

    // Ensure every filter in chrom_filters was used; if not, bail
    for (chrom, (_, _)) in filter_map.iter() {
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

    Ok(result)
}
