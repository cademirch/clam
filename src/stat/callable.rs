use anyhow::{Context, Ok, Result};
use d4::ssio::{D4TrackReader, D4TrackView};
use log::{debug, trace};
use std::collections::HashSet;
use std::fs::File;
use std::path::Path;
use d4::find_tracks;
pub struct D4CallableSites {
    readers: Vec<D4TrackReader<File>>,
}

impl D4CallableSites {
    pub fn from_file<P: AsRef<Path>>(pops: Option<Vec<&str>>, d4_file: P) -> Result<Self> {
        let d4_file_path = d4_file.as_ref();
        let mut readers: Vec<D4TrackReader<File>> = Vec::new();

        if let Some(pops) = pops {
            readers.reserve(pops.len());
            let mut found: Vec<&str> = vec![]; // Store found sample names as String
            let mut tracker = |p: Option<&Path>| {
                if let Some(path) = p {
                    let track_name = path.to_str().unwrap_or("");
                    // Check if any sample matches the track name
                    let this_pops = pops.clone();
                    if let Some(sample) = this_pops
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

            let file = File::open(&d4_file_path).unwrap();
            let mut found_tracks = vec![];
            find_tracks(file, &mut tracker, &mut found_tracks)?;
            for found in found_tracks {
                let file = File::open(&d4_file_path).unwrap();
                let rdr = D4TrackReader::from_reader(file, Some(&found.to_string_lossy())).unwrap();
                readers.push(rdr);
            }
        } else {
            let file = File::open(&d4_file_path)
                .context(format!("Couldn't open d4 file: {}", d4_file_path.display()))?;
            let rdr = D4TrackReader::from_reader(file, None)?;
            readers.push(rdr);
        }
        trace!("Created D4CallableSites");
        Ok(D4CallableSites { readers })
    }
    // pub fn query_counts(&mut self, chrom: &str, begin: u32, end: u32, ploidy: u32,) -> Result<Vec<Vec<u32>>> {
    //     // Initialize a 2D vector with default values of 0
    //     let adjusted_begin = begin - 1;
    //     let mut within_comps = vec![vec![0; (end - adjusted_begin) as usize]; self.readers.len()];
    
    //     // Create views for the specified region
    //     let mut views: Vec<D4TrackView<File>> = self
    //         .readers
    //         .iter_mut()
    //         .map(|reader| reader.get_view(chrom, adjusted_begin, end).unwrap())
    //         .collect();
    
    //     // Iterate through the range of positions
    //     for pos in adjusted_begin..end {
    //         // Get the callable individuals (values) for each population at this position
    //         let values: Vec<u32> = views
    //             .iter_mut()
    //             .map(|view| {
    //                 // Get the next reported position and value
    //                 let (reported_pos, value) = view.next().unwrap().unwrap();
    
    //                 // Ensure the reported position matches the current position
    //                 debug!(
    //                     "reported_pos: {}, pos: {}, value: {}",
    //                     reported_pos, pos, value
    //                 );
    //                 assert_eq!(reported_pos, pos);
    
    //                 value as u32
    //             })
    //             .collect();
    
    //         // Populate the results vector with values for each population
    //         for (pop_idx, &callable_indvs) in values.iter().enumerate() {
    //             let haps = ploidy * callable_indvs;
    //             let within_comp = if haps > 1 {
    //                 (haps as u64 * (haps as u64 - 1) / 2) as u32
    //             } else {
    //                 0
    //             };
    //             within_comps[pop_idx] += within_comp;
    //         }
    //     }
    
    //     Ok(res)
    // }
    
    pub fn query(
        &mut self,
        chrom: &str,
        begin: u32,
        end: u32,
        ploidy: u32,
        skip_sites: &HashSet<u32>,
    ) -> Result<(Vec<u32>, Option<Vec<u32>>)> {
        let num_pops = self.readers.len();
        let mut views: Vec<D4TrackView<File>> = self
            .readers
            .iter_mut()
            .map(|x| x.get_view(chrom, begin, end).unwrap())
            .collect();

        let mut within_comps: Vec<u32> = vec![0; num_pops];

        let num_pop_combs = if num_pops > 1 {
            num_pops * (num_pops - 1) / 2
        } else {
            0
        };

        let mut dxy_comps = if num_pop_combs > 0 {
            Some(vec![0; num_pop_combs])
        } else {
            None
        };

        for pos in begin..end {
            // Get the callable individuals (values) for each population at this position
            let values: Vec<u32> = views
                .iter_mut()
                .map(|view| {
                    let (reported_pos, value) = view.next().unwrap().unwrap();
                    
                    assert_eq!(reported_pos, pos);
                    value as u32
                })
                .collect();
            if skip_sites.contains(&pos) {
                // do this here so that view iterators are advanced still
                continue;
            }

            // Calculate within-population comparisons
            for (pop_idx, &callable_indvs) in values.iter().enumerate() {
                let haps = ploidy * callable_indvs;
                trace!("Ploidy: {}, Callable Indvs: {}, Callable haplotypes: {}", ploidy, callable_indvs, haps);
                let within_comp = if haps > 1 {
                    (haps as u64 * (haps as u64 - 1) / 2) as u32
                } else {
                    0
                };
                within_comps[pop_idx] += within_comp;
            }

            // Calculate between-population (dxy) comparisons
            if let Some(ref mut dxy_comps) = dxy_comps {
                let mut comb_idx = 0;
                for i in 0..num_pops {
                    let haps_i = ploidy * values[i];
                    for j in (i + 1)..num_pops {
                        let haps_j = ploidy * values[j];
                        dxy_comps[comb_idx] += haps_i * haps_j; // Pairwise comparison
                        comb_idx += 1;
                    }
                }
            }
        }
        debug!(
            "Callable Sites - {}:{}-{} {:?}",
            chrom, begin, end, within_comps
        );
        Ok((within_comps, dxy_comps))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use serde::Deserialize;
    use std::collections::HashSet;
    use std::fs::File;
    use std::path::Path;

    #[derive(Debug, Deserialize)]
    struct PopsTruthRecord {
        chrom: String,
        start: u32,
        end: u32,
        within_pop0: u32,
        within_pop1: u32,
        dxy: u32,
    }

    #[derive(Debug, Deserialize)]
    struct NoPopsTruthRecord {
        chrom: String,
        start: u32,
        end: u32,
        within_comp: u32,
    }

    fn init() {
        let _ = env_logger::builder()
            .target(env_logger::Target::Stdout)
            .filter_level(log::LevelFilter::Trace)
            .is_test(true)
            .try_init();
    }

    #[test]
    fn test_one_site() -> Result<()> {
        init();
        let fp = Path::new("tests/data/stat/diploid/no_pops_callable_sites.d4");
        let mut d4callable = D4CallableSites::from_file(None, fp)?;

        let (within_comps, _) = d4callable.query(
            "chr1",
            74,
            75,
            2,
            &HashSet::new(), // skip_sites, empty for this test
        )?;

        assert_eq!(within_comps[0], 45, "{:?}", within_comps);

        Ok(())
    }

    #[test]
    fn test_no_pops() -> Result<()> {
        init();
        let truth_fp = File::open("tests/data/stat/diploid/no_pops_truth.tsv")?;
        let mut reader = csv::ReaderBuilder::new()
            .delimiter(b'\t')
            .from_reader(truth_fp);
        let fp = Path::new("tests/data/stat/diploid/no_pops_callable_sites.d4");
        let mut d4callable = D4CallableSites::from_file(None, fp)?;

        for result in reader.deserialize() {
            let record: NoPopsTruthRecord = result?;
            let (within_comps, _) = d4callable.query(
                &record.chrom,
                record.start,
                record.end,
                2,
                &HashSet::new(), // skip_sites, empty for this test
            )?;

            assert_eq!(
                within_comps[0], record.within_comp,
                "Mismatch in within_comp for range {}:{}-{}",
                record.chrom, record.start, record.end
            );
        }
        Ok(())
    }

    #[test]
    fn test_no_pops_skip_sites() -> Result<()> {
        init();
        let truth_fp = File::open("tests/data/stat/diploid/no_pops_truth.tsv")?;
        let mut reader = csv::ReaderBuilder::new()
            .delimiter(b'\t')
            .from_reader(truth_fp);
        let fp = Path::new("tests/data/stat/diploid/no_pops_callable_sites.d4");
        let mut d4callable = D4CallableSites::from_file(None, fp)?;
        let mut skip_sites = HashSet::new();
        skip_sites.insert(10 as u32);
        for result in reader.deserialize() {
            let mut record: NoPopsTruthRecord = result?;
            let (within_comps, _) =
                d4callable.query(&record.chrom, record.start, record.end, 2, &skip_sites)?;
            if record.start <= 10 && record.end >= 10 {
                record.within_comp -= 120;
            }
            assert_eq!(
                within_comps[0], record.within_comp,
                "Mismatch in within_comp for range {}:{}-{}",
                record.chrom, record.start, record.end
            );
        }
        Ok(())
    }

    #[test]
    fn test_pops() -> Result<()> {
        init();
        let truth_fp = File::open("tests/data/stat/diploid/pops_truth.tsv")?;
        let mut reader = csv::ReaderBuilder::new()
            .delimiter(b'\t')
            .from_reader(truth_fp);
        let fp = Path::new("tests/data/stat/diploid/pops_callable_sites.d4");
        let pops = ["pop0", "pop1"].to_vec();
        let mut d4callable = D4CallableSites::from_file(Some(pops), fp)?;

        for result in reader.deserialize() {
            let record: PopsTruthRecord = result?;
            let (within_comps, dxy_comps) = d4callable.query(
                &record.chrom,
                record.start,
                record.end,
                2,
                &HashSet::new(), // skip_sites, empty for this test
            )?;

            assert_eq!(
                within_comps[0], record.within_pop0,
                "Mismatch in within_pop0 for range {}:{}-{}",
                record.chrom, record.start, record.end
            );
            assert_eq!(
                within_comps[1], record.within_pop1,
                "Mismatch in within_pop1 for range {}:{}-{}",
                record.chrom, record.start, record.end
            );

            assert!(dxy_comps.is_some());
            let dxy_comps = dxy_comps.unwrap();
            assert_eq!(
                dxy_comps[0], record.dxy,
                "Mismatch in dxy for range {}:{}-{}",
                record.chrom, record.start, record.end
            );
        }
        Ok(())
    }
}
