use std::panic;
use std::path::Path;
use std::time::Instant;

use anyhow::{bail, Context, Ok, Result};
use d4::task::{Task, TaskPartition};
use d4::{D4MatrixReader, D4TrackReader, MultiTrackReader};
use log::{info, trace};

use super::regions::CallableRegion;

pub struct TaskPart {
    parent_region: (String, u32, u32),
    thresholds: (f64, f64),
    count_sum: Vec<(u32, u32)>,
}

impl<T: Iterator<Item = i32> + ExactSizeIterator> TaskPartition<T> for TaskPart {
    type ParentType = TaskParent;
    type ResultType = Vec<(u32, u32)>;

    fn new(part_begin: u32, part_end: u32, parent: &Self::ParentType) -> Self {
        trace!(
            "Creating new TaskPart for region: {}-{}, in parent: {}:{}-{}.\nCountSum length: {}",
            part_begin,
            part_end,
            parent.chrom,
            parent.begin,
            parent.end,
            part_end - part_begin
        );
        Self {
            parent_region: (parent.chrom.clone(), parent.begin, parent.end),
            thresholds: parent.per_sample_thresholds,
            count_sum: vec![(0, 0); (part_end - part_begin) as usize],
        }
    }

    fn feed_range(&mut self, left: u32, right: u32, value: &mut T) -> bool {
        trace!(
            "TaskPart(parent: {:?}) Feeding range: {}-{}",
            self.parent_region,
            left,
            right
        );

        if right - left > 1 {
            panic!("Invalid range: {}-{} in task partition in parent: {:?}. The difference between right and left must be 1.", left, right, self.parent_region);
        }

        for v in value.into_iter() {
            trace!("Value: {:?}", v);
            let adjusted_left = (left as usize) % self.count_sum.len();

            self.count_sum[adjusted_left].1 += v as u32;
            trace!(
                "Testing if {:?} is greater than {} and less than {}",
                v,
                self.thresholds.0,
                self.thresholds.1
            );
            if v as f64 >= self.thresholds.0 && v as f64 <= self.thresholds.1 {
                self.count_sum[adjusted_left].0 += 1
            }
        }
        true
    }

    fn result(&mut self) -> Self::ResultType {
        trace!("Count sum: {:?}", self.count_sum);
        std::mem::take(&mut self.count_sum)
    }
}

pub struct TaskParent {
    chrom: String,
    begin: u32,
    end: u32,
    per_sample_thresholds: (f64, f64),
    mean_thresholds: (f64, f64),
    depth_proportion: f64,
    num_tracks: usize,
    output_counts: bool,
}

impl<T: Iterator<Item = i32> + ExactSizeIterator> Task<T> for TaskParent {
    type Partition = TaskPart;
    type Output = Vec<CallableRegion>;

    fn region(&self) -> (&str, u32, u32) {
        (self.chrom.as_str(), self.begin, self.end)
    }

    fn combine(&self, parts: &[Vec<(u32, u32)>]) -> Vec<CallableRegion> {
        trace!("Combining results from TaskParts: {:?}", parts);
        let mut res = Vec::new();
        let mut current_region: Option<CallableRegion> = None;

        for (index, (count, sum)) in parts.iter().flat_map(|v| v.iter()).enumerate() {
            trace!("count: {}", count);
            let mean = *sum as f64 / self.num_tracks as f64; // mean calculated based on number of individuals total
            let prop = *count as f64 / self.num_tracks as f64; // prop calculated based on number of individuals total
            trace!(
                "mean: {}, mean thresholds: {},{}",
                mean,
                self.mean_thresholds.0,
                self.mean_thresholds.1
            );
            if self.output_counts
                || ((self.mean_thresholds.0 < mean && mean < self.mean_thresholds.1)
                    && prop >= self.depth_proportion)
            {
                if let Some(ref mut region) = current_region {
                    // check if current_region is not none
                    // we have a current region, lets see if we can extend it
                    if region.end == (self.begin as usize + index) as u32
                        && (!self.output_counts || region.count == *count)
                    {
                        region.end += 1;
                    } else {
                        // cant extend
                        res.push(region.clone());
                        *region = CallableRegion {
                            count: *count,
                            begin: (self.begin as usize + index) as u32,
                            end: (self.begin as usize + index + 1) as u32,
                        };
                    }
                } else {
                    // no current region lets make one
                    current_region = Some(CallableRegion {
                        count: *count,
                        begin: (self.begin as usize + index) as u32,
                        end: (self.begin as usize + index + 1) as u32,
                    });
                }
            } else if let Some(region) = current_region.take() {
                res.push(region);
            }
        }

        // Push the last region if it exists
        if let Some(region) = current_region {
            res.push(region);
        }

        trace!("Combined result: {:?}", res);
        res
    }
}

pub fn prepare_tracks_from_file<P: AsRef<Path>>(
    d4_file_path: P,
    samples: Option<Vec<String>>,
) -> Result<Vec<D4TrackReader>> {
    let tracks = if let Some(samples) = samples {
        let mut found: Vec<String> = vec![]; // Store found sample names as String
        let tracker = |p: Option<&Path>| {
            if let Some(path) = p {
                let track_name = path.to_str().unwrap_or("");
                // Check if any sample matches the track name
                if let Some(sample) = samples
                    .iter()
                    .find(|sample| track_name.contains(sample.as_str()))
                {
                    found.push(sample.clone()); // Push the sample name to found as a String
                    true
                } else {
                    false
                }
            } else {
                false
            }
        };

        let tracks: Vec<D4TrackReader> =
            D4TrackReader::open_tracks(d4_file_path, tracker).context("Failed to open D4 file")?;

        let missing_samples: Vec<_> = samples
            .into_iter()
            .filter(|sample| !found.contains(sample))
            .collect();

        // If any samples are missing, return an error
        if !missing_samples.is_empty() {
            bail!("Samples not found in D4 file: {:?}", missing_samples);
        }
        tracks
    } else {
        // Get all tracks if no samples supplied
        let tracks: Vec<D4TrackReader> =
            D4TrackReader::open_tracks(d4_file_path, |_| true).context("Failed to open D4 file")?;
        tracks
    };

    Ok(tracks)
}

pub fn get_chrom_regions<P: AsRef<Path>>(
    d4_file_path: P,) -> Result<Vec<(String, u32, u32)>> {
        let tracks: Vec<D4TrackReader> =
            D4TrackReader::open_tracks(d4_file_path, |_| true).context("Failed to open D4 file")?;
            let reader = D4MatrixReader::new(tracks).unwrap();
        let regions = reader.chrom_regions().iter().map(|(chr, begin, end)| (chr.to_string(), *begin, *end)).collect();
        Ok(regions)
    }

pub fn run_tasks_on_tracks<P: AsRef<Path>>(
    d4_file_path: P,
    samples: Option<Vec<String>>,
    args: super::LociArgs
) -> Result<Vec<(String, u32, Vec<CallableRegion>)>> {
    let tracks = prepare_tracks_from_file(d4_file_path, samples)?;
    trace!("Successfully opened D4 file with {} tracks", tracks.len());
    let mean_thresholds = (args.mean_depth_min, args.mean_depth_max);
    let depth_proportion= args.depth_proportion;
    let output_counts = true;
    let num_tracks = tracks.len();
    let mut reader = D4MatrixReader::new(tracks).unwrap();

    let mut tasks = vec![];
    let regions = super::regions::prepare_chrom_regions(reader.chrom_regions(), args.get_per_sample_thresholds(), args.exclude_chrs.as_ref(), args.include_chrs.as_ref())?;
    for region in regions {
        tasks.push(TaskParent {
            chrom: region.chr.clone(),
            begin: region.begin,
            end: region.end,
            per_sample_thresholds: (region.min_filter, region.max_filter),
            mean_thresholds,
            depth_proportion,
            num_tracks,
            output_counts,
        });
    }

    let start = Instant::now();
    let result = reader.run_tasks(tasks)?;
    let duration = start.elapsed();
    info!("Ran D4 tasks in {:?}", duration);
    let formatted_result: Vec<(String, u32, Vec<CallableRegion>)> = result
        .into_iter()
        .map(|task| {
            let chrom_name = task.chrom.to_string(); // Get the chromosome name from TaskOutput
            let callable_regions = task.output.to_owned(); // Get the CallableRegion Vec from TaskOutput
            (chrom_name, 0, callable_regions) // Use 0 as the placeholder start position
        })
        .collect();
    Ok(formatted_result)
}

#[cfg(test)]
mod tests {

    use super::*;

    fn init() {
        let _ = env_logger::builder()
            .target(env_logger::Target::Stdout)
            .filter_level(log::LevelFilter::Trace)
            .is_test(true)
            .try_init();
    }

    type RowType = std::vec::IntoIter<i32>;
    fn create_task_parent(output_counts: bool) -> TaskParent {
        TaskParent {
            chrom: String::from("chr1"),
            begin: 0,
            end: 10,
            num_tracks: 10,
            mean_thresholds: (10.0, 100.0),
            depth_proportion: 0.5,
            per_sample_thresholds: (0.0, 100.0), // this doesnt matter for combine tests
            output_counts: output_counts,
        }
    }

    #[test]
    fn combine_output_counts_false_adjacent_regions_different_counts() {
        let output_counts = false;
        let task: TaskParent = create_task_parent(output_counts);

        let parts = vec![vec![(5, 101), (5, 101), (10, 101), (10, 101)]];

        let expected = vec![CallableRegion {
            count: 5,
            begin: 0,
            end: 4,
        }];

        let result: Vec<CallableRegion> = Task::<RowType>::combine(&task, &parts);

        assert_eq!(result, expected);
    }

    #[test]
    fn combine_output_counts_true_adjacent_regions_different_counts() {
        let output_counts = true;
        let task: TaskParent = create_task_parent(output_counts);

        let parts = vec![vec![(5, 101), (5, 101), (10, 101), (10, 101)]];

        let expected = vec![
            CallableRegion {
                count: 5,
                begin: 0,
                end: 2,
            },
            CallableRegion {
                count: 10,
                begin: 2,
                end: 4,
            },
        ];

        let result: Vec<CallableRegion> = Task::<RowType>::combine(&task, &parts);

        assert_eq!(result, expected);
    }

    #[test]
    fn combine_good_mean_bad_prop() {
        let output_counts = false;
        let task: TaskParent = create_task_parent(output_counts);

        let parts = vec![vec![(2, 101), (2, 101), (10, 101), (10, 101)]];

        let expected = vec![CallableRegion {
            count: 10,
            begin: 2,
            end: 4,
        }];

        let result: Vec<CallableRegion> = Task::<RowType>::combine(&task, &parts);

        assert_eq!(result, expected);
    }

    #[test]
    fn combine_single_continuous_region() {
        let output_counts = false;
        let task: TaskParent = create_task_parent(output_counts);

        let parts = vec![
            vec![(5, 50), (5, 50), (5, 50)], // mean depth is 5 not callable
        ];

        let expected = vec![];

        let result: Vec<CallableRegion> = Task::<RowType>::combine(&task, &parts);

        assert_eq!(result, expected);
    }

    #[test]
    fn combine_non_callable_regions() {
        let output_counts = false;
        let task: TaskParent = create_task_parent(output_counts);

        let parts = vec![
            vec![(1, 100), (1, 100), (5, 50), (5, 50)], // 0-2 okay mean but bad proportion. 2-4 okay proportion but bad mean
        ];

        let expected = vec![]; // No region should be returned

        let result: Vec<CallableRegion> = Task::<RowType>::combine(&task, &parts);

        assert_eq!(result, expected);
    }

    #[test]
    fn combine_multiple_separated_regions() {
        let output_counts = false;
        let task: TaskParent = create_task_parent(output_counts);

        let parts = vec![
            vec![(10, 101), (10, 101)], // callable
            vec![(5, 50), (5, 50)],     // not callable
            vec![(0, 0), (0, 0)],       // not callable
            vec![(10, 101), (10, 101)], // callable
        ];

        let expected = vec![
            CallableRegion {
                count: 10,
                begin: 0,
                end: 2,
            },
            CallableRegion {
                count: 10,
                begin: 6,
                end: 8,
            },
        ];

        let result: Vec<CallableRegion> = Task::<RowType>::combine(&task, &parts);

        assert_eq!(result, expected);
    }
    #[test]
    fn combine_below_mean_threshold() {
        let output_counts = false;
        let task: TaskParent = create_task_parent(output_counts);

        // Mean depth is 9.9, which is just below the lower threshold (10.0)
        let parts = vec![vec![(10, 99), (10, 99), (10, 99), (10, 99)]];

        let expected = vec![]; // No region should be returned since the mean is below the threshold

        let result: Vec<CallableRegion> = Task::<RowType>::combine(&task, &parts);

        assert_eq!(result, expected);
    }

    #[test]
    fn combine_above_mean_threshold() {
        init();
        let output_counts = false;
        let task: TaskParent = create_task_parent(output_counts);

        // Mean depth is 101, which is above the upper threshold (100.0)
        let parts = vec![vec![(10, 1001), (10, 1001), (10, 1001), (10, 1001)]];

        let expected = vec![]; // No region should be returned since the mean is above the threshold

        let result: Vec<CallableRegion> = Task::<RowType>::combine(&task, &parts);

        assert_eq!(result, expected);
    }

    #[test]
    fn combine_within_mean_threshold() {
        let output_counts = false;
        let task: TaskParent = create_task_parent(output_counts);

        // Mean depth is exactly 10, which is on the lower bound of the threshold
        let parts = vec![vec![(10, 101), (10, 101), (10, 999), (10, 999)]];

        let expected = vec![CallableRegion {
            count: 10,
            begin: 0,
            end: 4,
        }]; // The entire region should be returned since the mean is within the threshold

        let result: Vec<CallableRegion> = Task::<RowType>::combine(&task, &parts);

        assert_eq!(result, expected);
    }

    #[test]
    fn combine_on_upper_mean_threshold() {
        let output_counts = false;
        let task: TaskParent = create_task_parent(output_counts);

        let parts = vec![vec![(10, 999), (10, 999), (10, 999), (10, 999)]];

        let expected = vec![CallableRegion {
            count: 10,
            begin: 0,
            end: 4,
        }]; // The entire region should be returned since the mean is on the upper bound

        let result: Vec<CallableRegion> = Task::<RowType>::combine(&task, &parts);

        assert_eq!(result, expected);
    }

    #[test]
    fn combine_with_mean_exactly_on_boundaries() {
        init();
        let output_counts = false;
        let task: TaskParent = create_task_parent(output_counts);

        // Mean depth is exactly on the boundaries
        let parts = vec![
            vec![(10, 100), (10, 100)],   // mean = 10 (lower boundary)
            vec![(10, 1000), (10, 1000)], // mean = 100 (upper boundary)
        ];

        let expected = vec![];

        let result: Vec<CallableRegion> = Task::<RowType>::combine(&task, &parts);

        assert_eq!(result, expected);
    }

    #[test]
    fn combine_mean_below_threshold_non_callable() {
        let output_counts = false;
        let task: TaskParent = create_task_parent(output_counts);

        // Mean depth is too low, below the threshold
        let parts = vec![vec![(2, 5), (2, 5), (2, 5), (2, 5)]]; // mean depth = 2.5

        let expected = vec![]; // No region should be returned because the mean is too low

        let result: Vec<CallableRegion> = Task::<RowType>::combine(&task, &parts);

        assert_eq!(result, expected);
    }
}
