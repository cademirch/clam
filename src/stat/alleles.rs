use super::*;
use anyhow::{bail, Context, Result};
use crossbeam::channel::{bounded, Receiver, Sender};
use fnv::FnvHashMap;
use noodles::bgzf;
use noodles::core::Region;
use noodles::tabix;
use noodles::vcf::{
    self,
    variant::record::{
        samples::{keys::key, series::Value, Series},
        AlternateBases,
    },
};

use ::log::{debug, info};
use indicatif::{ProgressBar, ProgressStyle};
use std::num::NonZeroUsize;
use std::{fs::File, path::Path, thread};


pub type ProcessedRecord = (String, usize, usize, Option<Vec<[u32; 2]>>); //chrom idx, window_idx, position, allele counts

/// Count alleles in a single vcf record
pub fn count_alleles(
    record: &vcf::Record,
    header: &vcf::Header,
    sample_to_pop: &FnvHashMap<usize, usize>,
    counts: &mut Vec<[u32; 2]>,
) -> Result<()> {
    let samples = record.samples();
    let gt_series = samples
        .select(key::GENOTYPE)
        .ok_or_else(|| anyhow!("Malformed variant record: {:?}", record))?;

    for (sample_idx, value) in gt_series.iter(&header).enumerate() {
        let Some(Value::Genotype(genotype)) = value? else {
            bail!(
                "Invalid GT value for sample at index {}, in variant record: {:?}",
                sample_idx,
                record
            )
        };

        for result in genotype.iter() {
            let allele = result.with_context(|| {
                format!(
                    "Failed to parse alleles for sample at index {}, in variant record: {:?}",
                    sample_idx, record
                )
            })?;

            let (position, _) = allele;

            let pop_idx = sample_to_pop
                .get(&sample_idx)
                .ok_or_else(|| anyhow!("Unknown population for sample index {}", sample_idx))?;

            match position {
                Some(0) => counts[*pop_idx][0] += 1, // Reference allele
                Some(1) => counts[*pop_idx][1] += 1, // Alternate allele
                None => (),                          // Missing allele
                _ => bail!(
                    "Unexpected allele value in sample index {}: {:?}",
                    sample_idx,
                    position
                ),
            }
        }
    }

    Ok(())
}



// pub fn process_single_region<P: AsRef<Path>>(
//     region_idx: usize,
//     region: Region,
//     vcf: &P,
//     tbi: &P,
//     samp_idx_to_pop_idx: &FnvHashMap<usize, usize>,
//     processed_tx: &Sender<ProcessedRecord>,
// ) -> Result<()> {
//     let mut rdr = File::open(vcf.as_ref())
//         .map(bgzf::Reader::new)
//         .map(vcf::io::Reader::new)?;

//     let index = tabix::read(tbi.as_ref()).context("Failed to read .tbi file")?;

//     let header = rdr.read_header()?;

//     let query = rdr.query(&header, &index, &region)?;
//     trace!(
//         "Worker processing region {}, began iterating query",
//         region_idx
//     );
//     for result in query {
//         let record = result?;
//         let mut counts = None;

//         if record.alternate_bases().len() <= 1 {
//             let num_pops = samp_idx_to_pop_idx.values().max().map_or(0, |&v| v + 1);
//             let mut counts_vec: Vec<[u32; 2]> = vec![[0; 2]; num_pops];
//             count_alleles(&record, &header, &samp_idx_to_pop_idx, &mut counts_vec)?;
//             counts = Some(counts_vec);
//         }

//         let chrom = record.reference_sequence_name().to_string();
//         let start = record.variant_start().unwrap().unwrap().get();
//         processed_tx
//             .send((chrom, region_idx, start, counts))
//             .unwrap();
//     }

//     Ok(())
// }

// pub fn update_windows_for_record(
//     window_idx: usize,
//     start: usize,
//     counts: Option<Vec<[u32; 2]>>,
//     windows: &mut Vec<WindowedData>,
// ) -> Result<()> {
//     // Ensure the window index is within bounds
//     let windowed_data = windows
//         .get_mut(window_idx)
//         .ok_or_else(|| anyhow!("Window index {} is out of bounds", window_idx))?;

//     let should_update_counts = if let Some(sites) = &windowed_data.sites {
//         sites.contains(&(start as u32))
//     } else {
//         true // If `sites` is `None`, allow updating counts
//     };
//     trace!(
//         "Should update counts: {}, window idx: {}, startpos: {}, counts: {:?}",
//         should_update_counts,
//         window_idx,
//         start,
//         counts
//     );
//     if should_update_counts {
//         if let Some(counts) = counts {
//             // Update population counts
//             for (pop_idx, vals) in counts.iter().enumerate() {
//                 debug!("Updating popidx: {} with values: {:?}", pop_idx, vals);
//                 windowed_data.update_population(pop_idx, *vals);
//             }

//             // Update dxy differences if we have at least two populations
//             if counts.len() >= 2 {
//                 for pop1 in 0..counts.len() {
//                     for pop2 in (pop1 + 1)..counts.len() {
//                         let pop1_vals = counts[pop1];
//                         let pop2_vals = counts[pop2];

//                         let diffs = (pop1_vals[0] * pop2_vals[1]) + (pop1_vals[1] * pop2_vals[0]);
//                         let comps = (pop1_vals[0] + pop2_vals[1]) * (pop1_vals[1] + pop2_vals[0]);
//                         windowed_data.update_dxy_diffs(pop1, pop2, diffs, comps);
//                     }
//                 }
//             }
//         } else {
//             // If `counts` is `None`, mark `start` as skipped
//             windowed_data
//                 .sites_skipped
//                 .insert(start.try_into().unwrap());
//         }
//     }

//     Ok(())
// }
// pub fn process<P: AsRef<Path>>(
//     worker_count: NonZeroUsize,
//     regions: Vec<Region>,
//     vcf: P,
//     tbi: P,
//     samp_idx_to_pop_idx: FnvHashMap<usize, usize>,
//     mut windows: Vec<WindowedData>,
//     quiet: bool,
// ) -> Result<Vec<WindowedData>> {
//     let (region_tx, region_rx) = bounded(worker_count.get());
//     let (processed_tx, processed_rx) = bounded(worker_count.get());
//     let num_regions = regions.len();
//     info!("Processing {} regions.", num_regions);
//     let start_time = std::time::Instant::now();
//     let bar = if quiet {
//         None
//     } else {
//         let bar = ProgressBar::new(num_regions as u64);
//         bar.set_style(ProgressStyle::with_template("{spinner:.green} [{wide_bar:.cyan/blue}] {percent}% done.")
//         .unwrap()
//         .progress_chars("#>-"));
//     Some(bar)
    
//     };

//     thread::scope(|scope| {
//         scope.spawn(move || {
//             for (region_idx, region) in regions.into_iter().enumerate() {
//                 debug!("Sending region {} to region_tx", region_idx);
//                 region_tx.send((region_idx, region)).unwrap();
//             }
//             drop(region_tx);
//         });

//         for i in 0..worker_count.get() {
//             let region_rx = region_rx.clone();
//             let processed_tx = processed_tx.clone();
//             let vcf_path = vcf.as_ref().to_path_buf();
//             let tbi_path = tbi.as_ref().to_path_buf();
//             let samp_idx_to_pop_idx = samp_idx_to_pop_idx.clone();
//             let bar = bar.as_ref();

//             scope.spawn(move || -> Result<()> {
//                 let mut counter = 0;
//                 while let Ok((region_idx, region)) = region_rx.recv() {
//                     debug!("Worker {} processing region {}", i, region_idx);
//                     process_single_region(
//                         region_idx,
//                         region,
//                         &vcf_path,
//                         &tbi_path,
//                         &samp_idx_to_pop_idx,
//                         &processed_tx,
//                     )?;
//                     counter += 1;

//                     if counter >= num_regions / 100 {
//                         if let Some(bar) = bar {
//                             bar.inc((num_regions/100) as u64);
//                         }
//                         counter = 0;
//                     }
//                 }

//                 if let Some(bar) = bar {
//                     bar.inc(counter as u64);
//                 }

//                 Ok(())
//             });
//         }
//         drop(processed_tx);

//         while let Ok((_, window_idx, start, counts)) = processed_rx.recv() {
//             update_windows_for_record(window_idx, start, counts, &mut windows).unwrap();
//         }
//     });

//     if let Some(bar) = bar {
//         bar.finish();
//     }
//     info!("Processed {} regions in {:#?}", num_regions, start_time.elapsed());
//     Ok(windows)
// }

#[cfg(test)]
mod tests {
    use super::windows::*;
    use super::*;
    use anyhow::Result;
    use fnv::FnvHashMap;
    use noodles::core::Region;
    use noodles::vcf::io::Reader;
    use std::path::PathBuf;
    use std::str::FromStr;

    #[test]
    fn test_run() -> Result<()> {
        let vcf_path = PathBuf::from_str("tests/data/stat/diploid/small.vcf.gz");
        let tbi_path = PathBuf::from_str("tests/data/stat/diploid/small.vcf.gz.tbi");
        Ok(())
    }

    #[test]
    fn test_count_alleles_single_population() -> Result<()> {
        let data: &str = concat!(
            "##fileformat=VCFv4.3\n",
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample1\tsample2\n",
            "sq0\t1\t.\tA\tG\t.\tPASS\t.\tGT\t0/0\t0/1"
        );
        let data = data.as_bytes();

        let mut sample_to_pop = FnvHashMap::default();
        sample_to_pop.insert(0, 0);
        sample_to_pop.insert(1, 0);

        let mut reader = Reader::new(data);
        let actual_header = reader.read_header()?;

        for result in reader.records() {
            let record = result?;
            // Initialize counts with shape (num_pops, 5)
            let mut counts: Vec<[u32; 2]> = vec![[0; 2]; 1];

            // Pass the counts array to count_alleles
            count_alleles(&record, &actual_header, &sample_to_pop, &mut counts)?;

            assert_eq!(counts[0][0], 3); // ref alleles
            assert_eq!(counts[0][1], 1); // alt alleles
        }

        Ok(())
    }

    #[test]
    fn test_count_alleles_haploid_single_population() -> Result<()> {
        let data: &str = concat!(
            "##fileformat=VCFv4.3\n",
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample1\tsample2\n",
            "sq0\t1\t.\tA\tG\t.\tPASS\t.\tGT\t0\t1"
        );
        let data = data.as_bytes();

        let mut sample_to_pop = FnvHashMap::default();
        sample_to_pop.insert(0, 0);
        sample_to_pop.insert(1, 0);

        let mut reader = Reader::new(data);
        let actual_header = reader.read_header()?;

        for result in reader.records() {
            let record = result?;
            // Initialize counts with shape (num_pops, 5)
            let mut counts: Vec<[u32; 2]> = vec![[0; 2]; 1];

            // Pass the counts array to count_alleles
            count_alleles(&record, &actual_header, &sample_to_pop, &mut counts)?;

            assert_eq!(counts[0][0], 1); // ref alleles
            assert_eq!(counts[0][1], 1); // alt alleles
        }

        Ok(())
    }

    #[test]
    fn test_count_alleles_two_populations() -> Result<()> {
        let data: &str = concat!(
            "##fileformat=VCFv4.3\n",
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample1\tsample2\n",
            "sq0\t1\t.\tA\tG\t.\tPASS\t.\tGT\t0/0\t0/1"
        );
        let data = data.as_bytes();

        let mut sample_to_pop = FnvHashMap::default();
        sample_to_pop.insert(0, 0);
        sample_to_pop.insert(1, 1);

        let mut reader = Reader::new(data);
        let actual_header = reader.read_header()?;

        for result in reader.records() {
            let record = result?;
            // Initialize counts with shape (num_pops, 5)
            let mut counts: Vec<[u32; 2]> = vec![[0; 2]; 2];

            // Pass the counts array to count_alleles
            count_alleles(&record, &actual_header, &sample_to_pop, &mut counts)?;

            //pop0
            assert_eq!(counts[0][0], 2); // ref alleles
            assert_eq!(counts[0][1], 0); // alt alleles

            //pop1
            assert_eq!(counts[1][0], 1); // ref alleles
            assert_eq!(counts[1][1], 1); // alt alleles
        }

        Ok(())
    }
}
