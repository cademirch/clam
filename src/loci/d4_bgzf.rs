use super::{CallableRegion, ChromRegion};
use anyhow::{bail, Context, Result};
use crossbeam::channel::{bounded, Receiver, Sender};
use d4::{find_tracks, ssio::D4TrackReader, Chrom};
use log::trace;
use noodles::bgzf::{self, IndexedReader};
use rayon::prelude::*;
use std::time::Instant;
use std::{collections::HashMap, fs::File, path::Path, thread};

const CHUNK_SIZE: u32 = 10_000_000;

pub fn run_bgzf_tasks<P: AsRef<Path>>(
    d4_file: P,
    threads: usize,
    samples: Option<Vec<String>>,
    regions: Vec<ChromRegion>,
    mean_thresholds: (f64, f64),
    depth_proportion: f64,
    output_counts: bool,
) -> Result<Vec<(String, u32, Vec<CallableRegion>)>> {
    let (region_sender, region_receiver): (
        Sender<(String, u32, u32, f64, f64)>,
        Receiver<(String, u32, u32, f64, f64)>,
    ) = bounded(threads);
    let (result_sender, result_receiver) = bounded(threads);

    let region_sender_thread = {
        let region_sender = region_sender.clone();

        thread::spawn(move || {
            let chunked_regions: Vec<(String, u32, u32, f64, f64)> = regions
                .into_iter()
                .flat_map(|region| {
                    (region.begin..region.end).step_by(CHUNK_SIZE as usize).map(
                        move |chunk_start| {
                            let chunk_end = (chunk_start + CHUNK_SIZE).min(region.end);
                            (
                                region.chr.clone(),
                                chunk_start,
                                chunk_end,
                                region.min_filter,
                                region.max_filter,
                            )
                        },
                    )
                })
                .collect();

            // Send each chunked region to the channel
            for region in chunked_regions {
                region_sender.send(region).expect("Failed to send region");
            }
            drop(region_sender); // Close the sender when done
        })
    };

    let workers: Vec<_> = (0..threads)
        .map(|_| {
            let region_receiver = region_receiver.clone();
            let result_sender = result_sender.clone();
            let d4_file = d4_file.as_ref().to_path_buf();
            let samples = samples.clone();
            thread::spawn(move || {
                let track_names = samples.as_ref().map(|sample_names| {
                    sample_names
                        .iter()
                        .map(|s| s.as_str())
                        .collect::<Vec<&str>>()
                });

                let mut matrix_rdr = BgzfD4MatrixReader::from_path(d4_file, track_names)
                    .expect("Failed to create BgzfD4MatrixReader");
                let tid = thread::current().id();

                // Process each region received from the region channel
                for (chrom, begin, end, per_sample_min, per_sample_max) in region_receiver.iter() {
                    trace!(
                        "TID: {:?} | Processing region: {}:{}-{}",
                        tid,
                        chrom,
                        begin,
                        end
                    );
                    let start_time = Instant::now();
                    let callable = matrix_rdr
                        .get_callable_regions(
                            &chrom,
                            begin,
                            end,
                            (per_sample_min, per_sample_max),
                            mean_thresholds,
                            depth_proportion,
                            output_counts,
                        )
                        .expect("Failed to get callable regions");

                    // Send the result back through the result channel
                    let duration = start_time.elapsed();
                    trace!(
                        "TID: {:?} | Processed region: {}:{}-{} in {:?}",
                        tid,
                        chrom,
                        begin,
                        end,
                        duration
                    );
                    result_sender
                        .send((chrom, callable))
                        .expect("Failed to send result");
                }
            })
        })
        .collect();

    // Drop the original senders so worker threads will end when the region_sender_thread completes
    drop(region_sender);
    drop(result_sender);

    let mut callable_map: HashMap<String, Vec<CallableRegion>> = HashMap::new();
    for (chrom, callable_regions) in result_receiver {
        callable_map
            .entry(chrom)
            .or_insert_with(Vec::new)
            .extend(callable_regions);
    }

    // Wait for the region sender thread to complete
    region_sender_thread
        .join()
        .expect("Region sender thread panicked");

    // Wait for all worker threads to complete
    for worker in workers {
        worker.join().expect("Worker thread panicked");
    }
    for regions in callable_map.values_mut() {
        regions.sort_by_key(|region| region.begin);

        let mut merged_regions = Vec::new();
        let mut current_region: Option<CallableRegion> = None;

        for region in regions.drain(..) {
            if let Some(mut r) = current_region.take() {
                // Check if the current region is contiguous with the next
                if r.end == region.begin && r.count == region.count {
                    // Extend the current region
                    r.end = region.end;
                    current_region = Some(r);
                } else {
                    // Push the previous region and start a new one
                    merged_regions.push(r);
                    current_region = Some(region);
                }
            } else {
                // Start with the first region
                current_region = Some(region);
            }
        }
        // Push the final region
        if let Some(r) = current_region {
            merged_regions.push(r);
        }

        // Replace the old regions with merged regions
        *regions = merged_regions;
    }

    // Convert the merged map into the desired Vec format for BED output
    let regions: Vec<(String, u32, Vec<CallableRegion>)> = callable_map
        .into_iter()
        .map(|(chrom, callable_regions)| (chrom, 0, callable_regions)) // Use 0 as `begin` for each chromosome
        .collect();

    Ok(regions)
}

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
        let viewtime = Instant::now();
        let mut views: Vec<_> = self
            .tracks
            .iter_mut()
            .map(|track| track.get_view(chrom, begin, end).unwrap())
            .collect::<Vec<_>>();
        trace!(
            "Got {}:{}-{} views in {:?}",
            chrom,
            begin,
            end,
            viewtime.elapsed()
        );

        let mut callable_regions = Vec::new();
        let mut current_region: Option<CallableRegion> = None;
        let valtime = Instant::now();
        let cached_values: Vec<Vec<u32>> = views
            .par_iter_mut()
            .map(|view| {
                let mut pos_values = Vec::with_capacity((end - begin) as usize);
                for pos in begin..end {
                    let (reported_pos, value) = view.next().unwrap().unwrap();
                    assert_eq!(reported_pos, pos);
                    pos_values.push(value as u32);
                }
                pos_values
            })
            .collect();

        trace!(
            "Got values at {}:{}-{} in {:?}",
            chrom,
            begin,
            end,
            valtime.elapsed()
        );
        for pos_index in 0..(end - begin) as usize {
            // Collect values for each view at the current position in parallel
            let pos = begin + pos_index as u32;
            let values: Vec<u32> = cached_values.iter().map(|v| v[pos_index]).collect();
            trace!("chrom: {}, pos: {}, vals: {:?}", chrom, pos, values);
            // Calculate mean of values at this position
            let mean = values.iter().map(|&v| v as f64).sum::<f64>() / values.len() as f64;
            // trace!("Got values at pos {}", pos);
            // Count the values that pass per_sample_thresholds
            let count = values
                .iter()
                .filter(|&&v| {
                    let v = v as f64;
                    v >= per_sample_thresholds.0 && v <= per_sample_thresholds.1
                })
                .count() as u32;

            let meets_site_thresholds =
                mean > per_site_thresholds.0 && mean < per_site_thresholds.1;
            let meets_proportion = (count as f64) / (values.len() as f64) >= min_proportion;

            if output_counts || (meets_site_thresholds && meets_proportion) {
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


