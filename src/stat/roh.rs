use color_eyre::{
    eyre::{bail, Context},
    Result,
};
use log::debug;

use noodles::core::Region;
use noodles::tabix::io::indexed_reader::Builder as TabixReaderBuilder;
use rust_lapper::{Interval, Lapper};
use std::collections::HashMap;
use std::path::Path;

use crate::core::population::PopulationMap;

pub struct RohIndex {
    lapper: Lapper<u32, usize>,  
    pop_map: PopulationMap,
    sample_names: Vec<String>,
}

impl RohIndex {
    
    pub fn from_tabix_query(
        bed_path: &Path,
        region: &Region,
        analysis_samples: &[String],
        pop_map: PopulationMap,
    ) -> Result<Self> {
        let mut reader = TabixReaderBuilder::default()
            .build_from_path(bed_path)
            .wrap_err_with(|| format!("Failed to open ROH BED file: {}", bed_path.display()))?;
        
        let sample_names: Vec<String> = analysis_samples.to_vec();
        
        // Build lookup from sample name to index in analysis_samples
        let sample_lookup: HashMap<&str, usize> = analysis_samples
            .iter()
            .enumerate()
            .map(|(idx, name)| (name.as_str(), idx))
            .collect();
        
        let mut intervals: Vec<Interval<u32, usize>> = Vec::new();
        
        // Query the region
        let query = reader.query(region)
            .wrap_err_with(|| format!("Failed to query region: {:?}", region))?;
        
        for result in query {
            let record = result.wrap_err("Failed to read ROH record")?;
            let line = std::str::from_utf8(record.as_ref().as_bytes())
                .wrap_err("Invalid UTF-8 in ROH record")?;
            
            let fields: Vec<&str> = line.split('\t').collect();
            
            if fields.len() < 4 {
                bail!("ROH record missing required fields: {}", line);
            }
            
            let start: u32 = fields[1].parse()
                .wrap_err_with(|| format!("Invalid start position: {}", fields[1]))?;
            let end: u32 = fields[2].parse()
                .wrap_err_with(|| format!("Invalid end position: {}", fields[2]))?;
            let sample_name = fields[3];
            
            // Skip ROH entries for samples not in analysis set
            let sample_idx = match sample_lookup.get(sample_name) {
                Some(&idx) => idx,
                None => {
                    debug!(
                        "ROH entry for sample '{}' skipped (not in analysis samples)",
                        sample_name
                    );
                    continue;
                }
            };
            
            intervals.push(Interval {
                start: start + 1, // BED is 0-based, convert to 1-based
                stop: end + 1,    // Make half-open interval
                val: sample_idx,
            });
        }
        
        Ok(Self {
            lapper: Lapper::new(intervals),
            pop_map,
            sample_names,
        })
    }
    
    /// Check if a specific sample is in ROH at a given position
    #[inline]
    pub fn is_in_roh(&self, pos: usize, sample_idx: usize) -> bool {
        self.lapper
            .find(pos as u32, (pos + 1) as u32)
            .any(|iv| iv.val == sample_idx)
    }

    pub fn samples_in_roh_at(&self, pos: usize) -> Vec<usize> {
        self.lapper
            .find(pos as u32, (pos + 1) as u32)
            .map(|iv| iv.val)
            .collect()
    }
    
    /// Get ROH sample counts per population at a given position
    pub fn roh_counts_per_pop(&self, pos: usize) -> Vec<u16> {
        let mut counts = vec![0u16; self.pop_map.num_populations()];
        
        let samples_in_roh = self.samples_in_roh_at(pos);
        
        for &sample_idx in &samples_in_roh {
            if let Some(sample_name) = self.sample_names.get(sample_idx) {
                if let Some((pop_idx, _)) = self.pop_map.lookup(sample_name) {
                    counts[pop_idx] += 1;
                } else {
                    debug!(
                        "ROH entry for sample '{}' ignored (not in population map)",
                        sample_name
                    );
                }
            }
        }
        
        counts
    }
}