use super::counts::*;
use super::regions::{CallableRegion, ChromRegion};
use super::LociArgs;
use crate::stat::build_vcf_reader;
use anyhow::{bail, Result};
use indicatif::ProgressBar;
use log::warn;
use log::{info, trace};
use noodles::bgzf::Reader;
use noodles::core::Region;
use noodles::vcf::io::IndexedReader;
use noodles::vcf::variant::record::info::field::Value::Integer as info_int;
use noodles::vcf::variant::record::samples::series::Value::Integer as samples_int;
use noodles::vcf::variant::record::samples::Series;
use noodles::vcf::Header;
use std::collections::VecDeque;
use std::fs::File;
use std::path::{Path, PathBuf};
use std::sync::{Arc, Mutex};
use std::thread;
use std::time::Instant;
pub struct GvcfReader {
    inner: IndexedReader<Reader<File>>,
    header: Header,
    src: PathBuf,
}

impl GvcfReader {
    pub fn from_path<P: AsRef<Path>>(src: P) -> Result<Self> {
        let src_path = src.as_ref().to_path_buf();
        let (inner, header) = build_vcf_reader(src)?;
        Ok(Self {
            inner,
            header,
            src: src_path,
        })
    }
    pub fn chrom_regions(&self) -> Vec<(&str, u32, u32)> {
        let chroms = self.header.contigs();
        let mut res = Vec::with_capacity(chroms.len());

        for (name, contig) in chroms.iter() {
            if let Some(length) = contig.length() {
                res.push((
                    name.as_str(),
                    0,                                     // start position
                    length.try_into().unwrap_or(u32::MAX), // convert usize to u32
                ));
            }
        }

        res
    }
    pub fn get_depth(
        &mut self,
        chrom: &str,
        begin: u32,
        end: u32,
        depth_threshold: (f64, f64), // min and max depth thresholds
    ) -> Result<ChromosomeCounts> {
        let mut counts = ChromosomeCounts::new(chrom.to_string(), begin, end);
        let region: Region = format!("{}:{}-{}", chrom, begin, end).parse()?;
        let query = self.inner.query(&self.header, &region)?;
        trace!("Querying GVCF region {}:{}-{}", chrom, begin, end);
        for result in query {
            let record = result?;
            let Some(Ok(startpos)) = record.variant_start() else {
                warn!(
                    "Invalid start position in file {}\nOffending record: {:?}",
                    self.src.display(),
                    record
                );
                continue;
            };

            let info = record.info();
            let endpos = match info.get(&self.header, "END") {
                Some(Ok(Some(info_int(end)))) => end,
                _ => (startpos.get() + 1) as i32,
            };

            let samples = record.samples();
            let Some(format) = samples.select("DP") else {
                bail!(
                    "Missing DP format field in file {}\nOffending record: {:?}",
                    self.src.display(),
                    record
                );
            };
            let Some(Some(Ok(samples_int(dp)))) = format.get(&self.header, 0) else {
                bail!(
                    "Missing DP format field in file {}\nOffending record: {:?}",
                    self.src.display(),
                    record
                );
            };

            let array_start = (startpos.get() - 1) as usize;
            let array_end = (endpos - 1) as usize;
            trace!(
                "Record at {}:{}-{} DP={}, END={}",
                chrom,
                array_start + 1,
                array_end + 1,
                dp,
                endpos
            );
            for idx in array_start..=array_end {
                if idx < counts.counts.len() {
                    counts.depth_sums[idx] = dp as u32;
                    if (dp as f64) >= depth_threshold.0 && (dp as f64) <= depth_threshold.1 {
                        counts.counts.set(idx, true);
                    }
                }
            }
        }

        Ok(counts)
    }
}

pub fn run_tasks(
    gvcfs: VecDeque<PathBuf>,
    samples: Option<&[String]>,
    args: LociArgs,
    progress_bar: Option<ProgressBar>,
) -> Result<Vec<(String, u32, Vec<CallableRegion>)>> {
    let threads = args.threads;
    let mean_thresholds = (args.mean_depth_min, args.mean_depth_max);
    let min_proportion = args.depth_proportion;

    let files = if let Some(sample_names) = samples {
        gvcfs
            .into_iter()
            .filter(|path| {
                let path_str = path.to_string_lossy();
                sample_names.iter().any(|sample| path_str.contains(sample))
            })
            .collect()
    } else {
        gvcfs
    };
    let num_files = files.len();
    let regions = {
        let reader = GvcfReader::from_path(files.front().unwrap())?;
        super::regions::prepare_chrom_regions(
            reader.chrom_regions(),
            args.get_per_sample_thresholds(),
            args.exclude_chrs.as_ref(),
        )?
    };
    let work_queue = Arc::new(Mutex::new(files));

    let global_counts = Arc::new(Mutex::new(GlobalCounts::new(&regions)));

    thread::scope(|scope| {
        for i in 0..threads {
            trace!("Spawned worker {}", i);
            let work_queue = Arc::clone(&work_queue);
            let regions = &regions;
            let accumulator = Arc::clone(&global_counts);
            let bar = progress_bar.clone();
            scope.spawn(move || -> Result<()> {
                trace!("Worker {} starting...", i);
                loop {
                    let next_file = {
                        let mut queue = work_queue.lock().unwrap();
                        queue.pop_front()
                    };

                    match next_file {
                        Some(file) => {
                            let mut rdr = GvcfReader::from_path(file)?;
                            for region in regions {
                                let res = rdr.get_depth(
                                    &region.chr,
                                    region.begin,
                                    region.end,
                                    (region.min_filter, region.max_filter),
                                )?;
                                {
                                    let mut accumulator = accumulator.lock().unwrap();
                                    accumulator.merge(res);
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

    Ok(global_counts.finalize(min_proportion, mean_thresholds, num_files))
}
