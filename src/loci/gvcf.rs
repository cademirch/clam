use super::counts::*;
use super::regions::CallableRegion;
use super::LociArgs;
use crate::stat::build_vcf_reader;
use anyhow::{anyhow, Context, Result};
use indicatif::ProgressBar;

use log::{debug, error, info, trace, warn};
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
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::{Arc, Mutex};
use std::thread;
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
        min_gq: u32
    ) -> Result<ChromosomeCounts> {
        let mut counts = ChromosomeCounts::new(chrom.to_string(), begin, end);
        let region: Region = format!("{}:{}-{}", chrom, begin + 1, end).parse()?;
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
                debug!(
                    "Missing DP format field in file {}\nOffending record: {:?}",
                    self.src.display(),
                    record
                );
                continue;
            };
            let Some(Some(Ok(samples_int(dp)))) = format.get(&self.header, 0) else {
                debug!(
                    "Missing DP format field in file {}\nOffending record: {:?}",
                    self.src.display(),
                    record
                );
                continue;
            };

            let Some(format_gq) = samples.select("GQ") else {
                debug!(
                    "Missing GQ format field in file {}\nOffending record: {:?}",
                    self.src.display(),
                    record
                );
                continue;
            };
            let Some(Some(Ok(samples_int(gq)))) = format_gq.get(&self.header, 0) else {
                debug!(
                    "Missing GQ format field in file {}\nOffending record: {:?}",
                    self.src.display(),
                    record
                );
                continue;
            };
            let array_start = (startpos.get() - 1) as usize;
            let array_end = endpos as usize;
            trace!(
                "Record at {}:{}-{} DP={}, END={}, GQ={}",
                chrom,
                array_start + 1,
                array_end + 1,
                dp,
                endpos,
                gq
            );
            for idx in array_start..array_end {
                if idx < counts.counts.len() {
                    counts.depth_sums[idx] = dp as u32;
                    if gq >= min_gq as i32
                        && (dp as f64) >= depth_threshold.0
                        && (dp as f64) <= depth_threshold.1
                    {
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
    args: &LociArgs,
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

    info!("Processing {} files", num_files);
    let bar = if let Some(bar) = progress_bar {
        bar.set_length(num_files as u64);
        Some(bar)
    } else {
        None
    };

    let regions = {
        let reader = GvcfReader::from_path(files.front().context("No files to process")?)?;
        super::regions::prepare_chrom_regions(
            reader.chrom_regions(),
            args.get_per_sample_thresholds(),
            args.exclude_chrs.as_ref(),
            args.include_chrs.as_ref(),
        )?
    };

    let work_queue = Arc::new(Mutex::new(files));
    let global_counts = Arc::new(Mutex::new(GlobalCounts::new(&regions)));
    let completed = Arc::new(AtomicUsize::new(0));
    thread::scope(|scope| {
        let mut handles = Vec::new();

        for i in 0..threads {
            trace!("Spawned worker {}", i);
            let work_queue = Arc::clone(&work_queue);
            let regions = &regions;
            let accumulator = Arc::clone(&global_counts);
            let bar = bar.clone();
            let completed = Arc::clone(&completed);
            let handle = scope.spawn(move || -> Result<()> {
                trace!("Worker {} starting...", i);
                loop {
                    let next_file = {
                        let mut queue = work_queue.lock().unwrap();
                        queue.pop_front()
                    };
                    trace!("Worker {} got file: {:?}", i, next_file);

                    match next_file {
                        Some(file) => {
                            let start_time = std::time::Instant::now();
                            let mut rdr = GvcfReader::from_path(&file)?;
                            for region in regions {
                                trace!(
                                    "Worker {} processing region {}:{}-{}",
                                    i,
                                    region.chr,
                                    region.begin,
                                    region.end
                                );
                                match rdr.get_depth(
                                    &region.chr,
                                    region.begin,
                                    region.end,
                                    (region.min_filter, region.max_filter),
                                    args.min_gq
                                ) {
                                    Ok(res) => {
                                        let mut accumulator = accumulator.lock().unwrap();
                                        accumulator.merge(res);
                                    }
                                    Err(e) => {
                                        error!(
                                            "Error in get_depth for region {}:{}-{}: {}",
                                            region.chr, region.begin, region.end, e
                                        );
                                        return Err(e.into());
                                    }
                                }
                            }
                            if let Some(ref bar) = bar {
                                bar.inc(1);
                            }
                            let completed_count = completed.fetch_add(1, Ordering::SeqCst);
                            let remaining = num_files - completed_count - 1;
                            info!(
                                "Completed {} in {:.2?}. {} files remaining",
                                file.display(),
                                start_time.elapsed(),
                                remaining
                            );
                        }
                        None => break,
                    }
                }
                Ok(())
            });
            handles.push(handle);
        }

        for handle in handles {
            if let Err(e) = handle
                .join()
                .map_err(|e| anyhow!("Worker paniced: {:?}", e))?
            {
                return Err(e);
            }
        }
        Ok(())
    })?;

    let global_counts = Arc::try_unwrap(global_counts)
        .map_err(|_| anyhow::anyhow!("Failed to unwrap global counts"))?
        .into_inner()
        .map_err(|_| anyhow::anyhow!("Failed to get inner value of global counts"))?;

    Ok(global_counts.finalize(min_proportion, mean_thresholds, num_files))
}

#[cfg(test)]
mod tests {
    use crate::loci::regions::ChromRegion;

    #[test]
    fn test_gvcf() {
        let chroms = vec![
            (ChromRegion {
                chr: "NW_017805931.1".to_string(),
                begin: 0,
                end: 105,
                min_filter: 10.0,
                max_filter: f64::MAX,
            }),
        ];
        let infile = "tests/data/loci/test_no_pops/gvcf/s1.g.vcf.gz";
        let mut gvcf_reader = super::GvcfReader::from_path(infile).unwrap();
        let counts = gvcf_reader
            .get_depth("NW_017805931.1", 0, 105, (10.0, f64::MAX), 0)
            .unwrap();
        let mut global_counts = super::GlobalCounts::new(&chroms);
        global_counts.merge(counts);
        let callable = global_counts.finalize(0.0, (0.0, f64::MAX), 1);
        for c in callable {
            for j in c.2 {
                println!("{:?}", j);
            }
        }
    }
}
