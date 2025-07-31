use super::regions::CallableRegion;
use super::{counts::*, LociArgs};
use anyhow::{Context, Ok, Result};
use d4::find_tracks;
use d4::index::D4IndexCollection;
use d4::ssio::{D4TrackReader, D4TrackView};
use indicatif::ProgressBar;
use log::{trace, warn};
use noodles::bgzf::{self, IndexedReader};
use std::collections::VecDeque;
use std::fs::File;
use std::path::{Path, PathBuf};
use std::sync::{Arc, Mutex};
use std::thread;
use std::time::Instant;

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

    pub fn get_depth(
        &mut self,
        chrom: &str,
        begin: u32,
        end: u32,
        depth_threshold: (f64, f64), // min and max depth thresholds
    ) -> Result<ChromosomeCounts> {
        let mut view = self.inner_mut().get_view(chrom, begin, end)?;
        let mut counts = ChromosomeCounts::new(chrom.to_string(), begin, end);

        for pos in begin..end {
            let (reported_pos, value) = view.next().unwrap().unwrap();
            assert_eq!(reported_pos, pos);
            let depth: u32 = value.try_into().context("Couldn't fit depth into u32!")?;
            let idx = (pos - begin) as usize;
            counts.depth_sums[idx] = depth;
            if (value as f64) >= depth_threshold.0 && (value as f64) <= depth_threshold.1 {
                let mut bit = counts.counts.get_mut(idx).unwrap();
                *bit = true;
            }
        }

        Ok(counts)
    }
    pub fn chrom_regions(&self) -> Vec<(&str, u32, u32)> {
        self.inner()
            .get_header()
            .chrom_list()
            .iter()
            .map(|x| (x.name.as_str(), 0, x.size as u32))
            .collect()
    }
    pub fn get_callable_regions(
        &mut self,
        chrom: &str,
        begin: u32,
        end: u32,
        depth_threshold: (f64, f64), // min and max depth thresholds
        output_counts: bool,
    ) -> Result<Vec<CallableRegion>> {
        let mut view = self.inner_mut().get_view(chrom, begin, end)?;
        let mut callable_regions = Vec::new();
        let current_region: Option<CallableRegion> = None;
        let mut current_start = None;
        let mut current_value = None;
        let mut add_region = |start: u32, end: u32, value: f64| {
            let meets_threshold = value >= depth_threshold.0 && value <= depth_threshold.1;
            if meets_threshold || output_counts {
                callable_regions.push(CallableRegion {
                    begin: start,
                    end,
                    count: if meets_threshold { 1 } else { 0 },
                });
            }
        };

        // Process each position in the view
        for pos in begin..end {
            let (reported_pos, value) = view.next().unwrap().unwrap();
            assert_eq!(reported_pos, pos);
            let value = value as f64;
            match (current_start, current_value) {
                (Some(start), Some(prev_value)) if value == prev_value => {
                    // Continue current region
                    continue;
                }
                (Some(start), Some(prev_value)) => {
                    // End current region and maybe start new one
                    add_region(start, pos, prev_value);
                    current_start = Some(pos);
                    current_value = Some(value);
                }
                (None, None) => {
                    // Start new region
                    current_start = Some(pos);
                    current_value = Some(value);
                }
                _ => {}
            }
        }

        // Handle last region if exists
        if let (Some(start), Some(value)) = (current_start, current_value) {
            add_region(start, end, value);
        }

        Ok(callable_regions)
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

pub fn run_tasks(
    work_queue: VecDeque<(PathBuf, Option<usize>)>,
    args: LociArgs,
    progress_bar: Option<ProgressBar>,
    regions: Vec<super::regions::ChromRegion>
) -> Result<(
    Vec<(String, u32, Vec<CallableRegion>)>,      // global
    Vec<Vec<(String, u32, Vec<CallableRegion>)>>, // per-population
)> {
    let threads = args.threads;
    let mean_thresholds = (args.mean_depth_min, args.mean_depth_max);
    let min_proportion = args.depth_proportion;
    let total_num_samples = work_queue.len();

    

    let global_counts = Arc::new(Mutex::new(GlobalCounts::new(&regions)));

    let (pop_accumulator, pop_map) = if let Some(pop_map) = &args.population_map {
    let pop_accumulator = Arc::new(Mutex::new(
        (0..pop_map.num_populations())
            .map(|_| GlobalCounts::new(&regions))
            .collect::<Vec<_>>(),
    ));
    (Some(pop_accumulator), Some(pop_map))
} else {
    (None, None)
};
    let work_queue = Arc::new(Mutex::new(work_queue));

    thread::scope(|scope| {
        for i in 0..threads {
            let work_queue = Arc::clone(&work_queue);
            let regions = &regions;
            let pop_accumulator = pop_accumulator.as_ref().map(Arc::clone);
            let accumulator = Arc::clone(&global_counts);
            let bar = progress_bar.clone();
            scope.spawn(move || -> Result<()> {
                loop {
                    let next_file = {
                        let mut queue = work_queue.lock().unwrap();
                        queue.pop_front()
                    };
                    
                    match next_file {
                        Some((file, pop_idx_opt)) => {
                            let mut rdr = BGZID4TrackReader::from_path(file, None)?;
                            for region in regions {
                                let res = rdr.get_depth(
                                    &region.chr,
                                    region.begin,
                                    region.end,
                                    (region.min_filter, region.max_filter),
                                )?;
                                
                                {
                                    let mut accumulator = accumulator.lock().unwrap();
                                    accumulator.merge(res.clone());
                                }
                                if let (Some(pop_acc), Some(pop_idx)) =
                                    (pop_accumulator.as_ref(), pop_idx_opt)
                                {
                                    let mut pop_acc = pop_acc.lock().unwrap();
                                    pop_acc[pop_idx].merge(res);
                                }
                            }
                            if let Some(ref bar) = bar {
                                bar.inc(1);
                            }
                        }
                        None => break,
                    }
                }
                Ok(())
            });
        }
    });

    let global_counts = Arc::try_unwrap(global_counts)
        .map_err(|_| anyhow::anyhow!("Failed to unwrap global counts"))?
        .into_inner()
        .map_err(|_| anyhow::anyhow!("Failed to get inner value of global counts"))?;

    let finalized_pop_counts =
        if let (Some(pop_accumulator), Some(pop_map)) = (pop_accumulator, pop_map) {
            let pop_counts = Arc::try_unwrap(pop_accumulator)
                .map_err(|_| anyhow::anyhow!("Failed to unwrap pop_accumulator"))?
                .into_inner()
                .map_err(|_| anyhow::anyhow!("Failed to get inner value of pop_accumulator"))?;
            let sample_counts = pop_map.get_sample_counts_per_population();
            pop_counts
                .into_iter()
                .zip(sample_counts.into_iter())
                .map(|(pop_gc, num_samples)| pop_gc.finalize(0.0, mean_thresholds, num_samples))
                .collect()
        } else {
            Vec::new()
        };

    Ok((
        global_counts.finalize(min_proportion, mean_thresholds, total_num_samples),
        finalized_pop_counts,
    ))
}

// pub fn run_bgzf_tasks<P: AsRef<Path>>(
//     d4_file: P,
//     threads: usize,
//     samples: Option<Vec<String>>,
//     regions: Vec<ChromRegion>,
//     depth_thresholds: (f64, f64),
//     min_proportion: f64,
//     output_counts: bool,
// ) -> Result<Vec<(String, u32, Vec<CallableRegion>)>> {
//     // Initialize global counts
//     let global_counts = Arc::new(Mutex::new(GlobalCounts::new(&regions)));

//     let d4_path = d4_file.as_ref();

//     // Get list of tracks/files to process
//     let file_queue: VecDeque<(PathBuf, Option<String>)> = if let Some(extension) = d4_path.extension().and_then(|e| e.to_str()) {
//         match extension {
//             "gz" => {
//                 // For merged files, we'll read tracks from the same file but with different track names
//                 let indexed_reader = build_bgzf_reader(d4_path)?;
//                 let mut track_paths = vec![];
//                 find_tracks(indexed_reader, |_| true, &mut track_paths)?;

//                 let tracks: VecDeque<_> = if let Some(ref sample_names) = samples {
//                     track_paths.into_iter()
//                         .filter_map(|track_path| {
//                             let track_str = track_path.to_string_lossy();
//                             if sample_names.iter().any(|sample| track_str.contains(sample)) {
//                                 Some((d4_path.to_path_buf(), Some(track_str.into_owned())))
//                             } else {
//                                 None
//                             }
//                         })
//                         .collect()
//                 } else {
//                     track_paths.into_iter()
//                         .map(|track_path| {
//                             (d4_path.to_path_buf(), Some(track_path.to_string_lossy().into_owned()))
//                         })
//                         .collect()
//                 };
//                 tracks
//             },
//             "txt" => {
//                 // For file of files, read separate files
//                 let paths: VecDeque<(PathBuf, Option<String>)> = std::fs::read_to_string(d4_path)?
//                     .lines()
//                     .map(|line| (PathBuf::from(line.trim()), None))
//                     .filter(|(path, _)| {
//                         if let Some(ref sample_names) = samples {
//                             let path_str = path.to_string_lossy();
//                             sample_names.iter().any(|sample| path_str.contains(sample))
//                         } else {
//                             true
//                         }
//                     })
//                     .collect();
//                 paths
//             },
//             _ => return Err(anyhow::anyhow!("Unsupported file extension")),
//         }
//     } else {
//         return Err(anyhow::anyhow!("File has no extension"));
//     };

//     let num_files = file_queue.len();
//     info!("Processing {} files", num_files);

//     // Create shared work queue
//     let work_queue = Arc::new(Mutex::new(file_queue));

//     // Track progress
//     let start_time = Instant::now();
//     let files_processed = Arc::new(Mutex::new(0usize));
//     let last_log_time = Arc::new(Mutex::new(Instant::now()));

//     // Use scope to ensure all threads complete before we finalize results
//     thread::scope(|scope| {
//         let mut handles = Vec::new();

//         // Spawn worker threads
//         for thread_id in 0..threads {
//             let work_queue = Arc::clone(&work_queue);
//             let global_counts = Arc::clone(&global_counts);
//             let files_processed = Arc::clone(&files_processed);
//             let last_log_time = Arc::clone(&last_log_time);
//             let regions = regions.clone();

//             handles.push(scope.spawn(move || {
//                 loop {
//                     // Get next file from queue
//                     let work_item = {
//                         let mut queue = work_queue.lock().unwrap();
//                         queue.pop_front()
//                     };

//                     // Break if no more files
//                     let (file_path, track_name) = match work_item {
//                         Some(item) => item,
//                         None => break,
//                     };

//                     trace!("Thread {} processing file: {} track: {:?}",
//                           thread_id, file_path.display(), track_name);

//                     // Process the file
//                     match process_single_file(&file_path, track_name.as_deref(), &regions, depth_thresholds) {
//                         Ok(sample_regions) => {
//                             // Add results to global counts
//                             let mut global = global_counts.lock().unwrap();
//                             // global.add_sample_regions(sample_regions);

//                             // Update progress
//                             let mut processed = files_processed.lock().unwrap();
//                             *processed += 1;

//                             // Log progress every 5 files or at completion
//                             if *processed % 5 == 0 || *processed == num_files {
//                                 let mut last_time = last_log_time.lock().unwrap();
//                                 let elapsed_since_last = last_time.elapsed();
//                                 *last_time = Instant::now();

//                                 info!(
//                                     "Processed {}/{} files ({:.1}%) | Last 5 files in {:?} | Total time: {:?}",
//                                     processed,
//                                     num_files,
//                                     (*processed as f64 / num_files as f64) * 100.0,
//                                     elapsed_since_last,
//                                     start_time.elapsed()
//                                 );
//                             }
//                         },
//                         Err(e) => {
//                             warn!("Failed to process file {}: {}", file_path.display(), e);
//                         }
//                     }
//                 }
//             }));
//         }

//         // Wait for all threads to complete
//         for handle in handles {
//             handle.join().unwrap();
//         }
//     });

//     // Finalize results
//     let global_counts = Arc::try_unwrap(global_counts)
//         .map_err(|_| anyhow::anyhow!("Failed to unwrap global counts"))?
//         .into_inner()
//         .map_err(|_| anyhow::anyhow!("Failed to get inner value of global counts"))?;

//     Ok(global_counts.finalize(min_proportion, 0.0, 100))
// }

// fn process_single_file(
//     file_path: &Path,
//     track_name: Option<&str>,
//     regions: &[ChromRegion],
//     depth_thresholds: (f64, f64),
// ) -> Result<Vec<(String, Vec<CallableRegion>)>> {
//     let mut reader = BGZID4TrackReader::from_path(file_path, track_name)?;
//     let mut results = Vec::new();

//     for region in regions {
//         let callable_regions = reader.get_callable_regions(
//             &region.chr,
//             region.begin,
//             region.end,
//             depth_thresholds,
//             false,  // output_counts handled at merge time
//         )?;

//         if !callable_regions.is_empty() {
//             results.push((region.chr.clone(), callable_regions));
//         }
//     }

//     Ok(results)
// }
