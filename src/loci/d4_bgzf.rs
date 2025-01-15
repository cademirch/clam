use std::collections::{HashSet,HashMap};
use std::fs::File;
use std::path::Path;
use std::thread;
use std::time::Instant;

use anyhow::{Context, Result};
use crossbeam::channel::{bounded, Receiver, Sender};
use d4::find_tracks;
use d4::index::D4IndexCollection;
use d4::ssio::{D4TrackReader, D4TrackView};
use log::{info, trace, warn};
use noodles::bgzf::{self, IndexedReader};

use super::{CallableRegion, ChromRegion};

const CHUNK_SIZE: u32 = 10_000_000;

/// A track reader for a single file/track.
/// File must be bgzipped with a `.gzi` index file next to it.
///
fn build_bgzf_reader<P: AsRef<Path>>(src: P) -> Result<IndexedReader<File>> {
    bgzf::indexed_reader::Builder::default()
        .build_from_path(src.as_ref())
        .with_context(|| format!("Failed to read d4 file: {}", src.as_ref().display()))
}
pub struct BGZID4TrackReader {
    inner: D4TrackReader<IndexedReader<File>>,
}

impl BGZID4TrackReader {
    pub fn from_path<P: AsRef<Path>>(src: P, track_name: Option<&str>) -> Result<Self> {
        let indexed_reader = build_bgzf_reader(src.as_ref())?;

        let d4_reader = D4TrackReader::from_reader(indexed_reader, track_name).map_err(|e| {
            anyhow::anyhow!(
                "Failed to create D4 track reader for file: {}, track_name: {:?}. Error: {}",
                src.as_ref().display(),
                track_name,
                e
            )
        })?;
        let track_root = d4_reader.as_root();
        let ic = D4IndexCollection::from_root_container(&track_root);

        if ic.is_err() {
            warn!("SFI index not found for D4 file: {}, this will result in slower performance. You can create the index by running d4tools index build {}", &src.as_ref().display(), &src.as_ref().display())
        }

        Ok(Self { inner: d4_reader })
    }
    pub fn inner_mut(&mut self) -> &mut D4TrackReader<IndexedReader<File>> {
        &mut self.inner
    }
    pub fn inner(&self) -> &D4TrackReader<IndexedReader<File>> {
        &self.inner
    }
}

pub struct BGZID4MatrixReader {
    readers: Vec<BGZID4TrackReader>,
}

impl BGZID4MatrixReader {
    pub fn from_paths<P: AsRef<Path>>(paths: Vec<P>) -> Result<Self> {
        // Create readers from the provided paths
        let readers: Result<Vec<_>> = paths
            .into_iter()
            .map(|p| BGZID4TrackReader::from_path(p, None))
            .collect();
    
        Ok(Self { readers: readers? })
    }
    pub fn from_merged<P: AsRef<Path>>(src: P, track_names: Option<Vec<&str>>) -> Result<Self> {
        let indexed_reader = build_bgzf_reader(src.as_ref())?;

        if let Some(track_names) = track_names {
            let mut found: Vec<&str> = vec![]; // Store found sample names as String
            let tracker = |p: Option<&Path>| {
                if let Some(path) = p {
                    let track_name = path.to_str().unwrap_or("");
                    // Check if any sample matches the track name
                    if let Some(sample) = track_names
                        .iter()
                        .find(|&&sample| track_name.contains(sample))
                    {
                        found.push(sample); // Push the sample name to found as a String
                        true
                    } else {
                        false
                    }
                } else {
                    false
                }
            };
            let mut found_tracks = vec![];

            find_tracks(indexed_reader, tracker, &mut found_tracks)?;
            let readers: Result<Vec<_>> = found_tracks
                .into_iter()
                .map(|track| {
                    BGZID4TrackReader::from_path(src.as_ref(), Some(&track.to_string_lossy()))
                })
                .collect();

            Ok(Self { readers: readers? })
        } else {
            let mut track_names = vec![];
            find_tracks(indexed_reader, |_| true, &mut track_names)?;
            let readers: Result<Vec<_>> = track_names
                .into_iter()
                .map(|track| {
                    BGZID4TrackReader::from_path(src.as_ref(), Some(&track.to_string_lossy()))
                })
                .collect();

            Ok(Self { readers: readers? })
        }
    }
    pub fn chrom_regions(&self) -> Vec<(&str, u32, u32)> {
        self.readers[0]
            .inner()
            .get_header()
            .chrom_list()
            .iter()
            .map(|x| (x.name.as_str(), 0, x.size as u32))
            .collect()
    }
    pub fn get_views(
        &mut self,
        chrom: &str,
        begin: u32,
        end: u32,
    ) -> Result<Vec<D4TrackView<IndexedReader<File>>>> {
        let views: Vec<_> = self
            .readers
            .iter_mut()
            .map(|track| track.inner_mut().get_view(chrom, begin, end).unwrap())
            .collect::<Vec<_>>();
        Ok(views)
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
        let viewtime = Instant::now();
        let mut views = self.get_views(chrom, begin, end)?;
        trace!("Got views in {:?}", viewtime.elapsed());
        let mut callable_regions = Vec::new();
        let mut current_region: Option<CallableRegion> = None;

        let cached_values: Vec<Vec<u32>> = views
            .iter_mut()
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
        for pos_index in 0..(end - begin) as usize {
            // Collect values for each view at the current position in parallel
            let pos = begin + pos_index as u32;
            let values: Vec<u32> = cached_values.iter().map(|v| v[pos_index]).collect();
            // trace!("chrom: {}, pos: {}, vals: {:?}", chrom, pos, values);
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
    let total_chunks: u32 = regions
        .iter()
        .map(|r| (r.end - r.begin + CHUNK_SIZE - 1) / CHUNK_SIZE) // Ceiling division to count all chunks
        .sum();
    info!("Total chunks to process: {}", total_chunks);
    let mut completed_chunks = 0;

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
    let start_time = Instant::now();
    let mut last_log_time = Instant::now();

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

                let mut matrix_rdr = if let Some(extension) = Path::new(&d4_file).extension().and_then(|e| e.to_str()) {
                    match extension {
                        "gz" => {
                            BGZID4MatrixReader::from_merged(d4_file, track_names)
                                .expect("Failed to create BgzfD4MatrixReader from merged")
                        }
                        // TODO write function to handle this and do error handling n stuff
                        "txt" => {
                            // Read paths from the .txt file
                            let paths = std::fs::read_to_string(&d4_file)
                                .expect("Failed to read paths from txt file")
                                .lines()
                                .map(|line| line.trim().to_string())
                                .collect::<Vec<_>>();
                            trace!("FOF Paths: {:?}", paths);
                            // Filter paths based on track_names if provided
                            let filtered_paths: Vec<_> = if let Some(track_names) = &track_names {
                                trace!("FOF Track Names: {:?}", track_names);
                                let track_name_set: HashSet<_> = track_names.iter().cloned().collect();
                                paths
                                    .into_iter()
                                    .filter(|path| {
                                        track_name_set.iter().any(|track_name| path.contains(track_name))
                                    })
                                    .collect()
                            } else {
                                paths
                            };
                            
                            trace!("Filtered paths: {:?}", filtered_paths);
                            BGZID4MatrixReader::from_paths(filtered_paths)
                                .expect("Failed to create BgzfD4MatrixReader from paths")
                        }
                        _ => panic!("Unsupported file extension: {}", extension),
                    }
                } else {
                    panic!("File does not have a valid extension");
                };
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
        completed_chunks += 1;

        // Log progress every 10 chunks
        if completed_chunks % 10 == 0 {
            let elapsed_for_last_10 = last_log_time.elapsed();
            last_log_time = Instant::now(); // Reset for the next 10 chunks
            let total_elapsed = start_time.elapsed();

            info!(
                "Processed {} chunks out of {} ({:.2}%) | Last 10 chunks in {:?} | Total time: {:?}",
                completed_chunks,
                total_chunks,
                completed_chunks as f64 / total_chunks as f64 * 100.0,
                elapsed_for_last_10,
                total_elapsed
            );
        }
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
