use anyhow::{bail, Context, Result};
use d4::{
    find_tracks,
    ssio::{D4TrackReader, D4TrackView},
    Chrom,
};
use noodles::bgzf::{self, IndexedReader};
use rayon::prelude::*;
use std::{
    collections::HashSet,
    fs::File,
    io::{BufReader, Read, Seek},
    path::Path,
    thread,
};
use crate::d4_tasks::CallableRegion;

pub struct BgzfD4MatrixReader {
    tracks: Vec<D4TrackReader<IndexedReader<File>>>,
}

impl BgzfD4MatrixReader {
    pub fn from_path<P: AsRef<Path>>(d4_file: P, samples: Option<Vec<&str>>) -> Result<Self> {
        let mut bgzf_reader = bgzf::indexed_reader::Builder::default()
            .build_from_path(d4_file.as_ref())
            .context("Failed to open d4 file!")?;
        let mut found = vec![];
        find_tracks(bgzf_reader, |_| true, &mut found)?; // just get all tracks

        if let Some(samples_vec) = &samples {
            let mut missing_samples = Vec::new();
            let mut filtered_found = Vec::new();

            for sample in samples_vec {
                let mut is_sample_found = false;
                for track_name in &found {
                    if track_name.to_string_lossy().contains(sample) {
                        filtered_found.push(track_name.clone());
                        is_sample_found = true;
                    }
                }
                if !is_sample_found {
                    missing_samples.push(sample.to_string());
                }
            }

            if !missing_samples.is_empty() {
                bail!("Samples not found in d4 file! {:?}", missing_samples);
            }

            found = filtered_found;
        }
        let mut readers = Vec::new();
        for track_name in &found {
            let bgzf_reader = bgzf::indexed_reader::Builder::default()
                .build_from_path(d4_file.as_ref())
                .context("Failed to open d4 file!")?;

            let rdr = D4TrackReader::from_reader(bgzf_reader, Some(&track_name.to_string_lossy()))
                .context(format!(
                    "Failed to create D4TrackReader for track {}",
                    track_name.display()
                ))?;
            readers.push(rdr);
        }

        Ok(Self { tracks: readers })
    }
    pub fn get_view(
        &mut self,
        chrom: &str,
        begin: u32,
        end: u32,
        per_sample_thresholds:(f64,f64),
        per_site_thresholds:(f64,f64),
        min_proportion: f64
    ) -> Result<Vec<Vec<u32>>> {
        let mut views: Vec<_> = self.tracks
        .iter_mut()
        .map(|track| track.get_view(chrom, begin, end).unwrap())
        .collect();

    // Parallelize over the positions in the specified range
    let mut all_values = Vec::with_capacity((end - begin) as usize);
    
    // Process each position sequentially
    for pos in begin..end {
        // Collect values for each view at the current position in parallel
        let values: Vec<u32> = views
            .par_iter_mut()
            .map(|view| {
                let (reported_pos, value) = view.next().unwrap().unwrap();
                assert_eq!(reported_pos, pos);
                value as u32
            })
            .collect();

        // Add the collected values for this position to the result
        all_values.push(values);
    }

    Ok(all_values)
    }
    pub fn chrom_regions(&self) -> Vec<(&str, u32, u32)> {
        self.tracks[0]
            .get_header()
            .chrom_list()
            .iter()
            .map(|x| (x.name.as_str(), 0, x.size as u32))
            .collect()
    }
    pub fn chrom_list(&self) -> Vec<&Chrom> {
        self.tracks[0].get_header().chrom_list().iter().collect()
    }
    pub fn get_callable_regions(
        &mut self,
        chrom: &str,
        begin: u32,
        end: u32,
        per_sample_thresholds: (f64, f64),
        per_site_thresholds: (f64, f64),
        min_proportion: f64,
        output_counts: bool,
    ) -> Result<Vec<CallableRegion>> {
        // Prepare views for processing
        let mut views: Vec<_> = self
            .tracks
            .iter_mut()
            .map(|track| track.get_view(chrom, begin, end).unwrap()) 
            .collect::<Vec<_>>();
    
        
        let mut callable_regions = Vec::new();
        let mut current_region: Option<CallableRegion> = None;
    
        
        for pos in begin..end {
            // Collect values for each view at the current position in parallel
            let values: Vec<u32> = views
                .par_iter_mut()
                .map(|view| {
                    let (reported_pos, value) = view.next().unwrap().unwrap();
                    assert_eq!(reported_pos, pos);
                    value as u32
                })
                .collect();
    
            // Calculate mean of values at this position
            let mean = values.iter().map(|&v| v as f64).sum::<f64>() / values.len() as f64;
    
            // Count the values that pass per_sample_thresholds
            let count = values
                .iter()
                .filter(|&&v| {
                    let v = v as f64;
                    v >= per_sample_thresholds.0 && v <= per_sample_thresholds.1
                })
                .count() as u32;
    
            
            let meets_site_thresholds = mean >= per_site_thresholds.0 && mean <= per_site_thresholds.1;
            let meets_proportion = (count as f64) / (values.len() as f64) >= min_proportion;
    
            if meets_site_thresholds && meets_proportion {
                // Start a new region or extend the current one
                if let Some(ref mut region) = current_region {
                    // Check if we can extend the region
                    if region.end == pos && (!output_counts || region.count == count) {
                        region.end += 1; // Extend the current region
                    } else {
                        // Can't extend, finalize the region and start a new one
                        callable_regions.push(region.clone());
                        *region = CallableRegion {
                            count,
                            begin: pos,
                            end: pos + 1,
                        };
                    }
                } else {
                    // No current region, so start a new one
                    current_region = Some(CallableRegion {
                        count,
                        begin: pos,
                        end: pos + 1,
                    });
                }
            } else if let Some(region) = current_region.take() {
                // Finalize any open region if criteria are no longer met
                callable_regions.push(region);
            }
        }
    
        // If there's an ongoing region at the end of the loop, finalize it
        if let Some(region) = current_region {
            callable_regions.push(region);
        }
    
        Ok(callable_regions)
    }
}
