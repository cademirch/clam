use std::collections::HashSet;
use std::fs::File;
use std::path::Path;

use anyhow::{Context, Ok, Result};
use bstr::BString;
use d4::find_tracks;
use d4::ssio::{D4TrackReader, D4TrackView};
use indexmap::IndexSet;
use log::{debug, trace};
use noodles::bed;
use noodles::bgzf::{self as bgzf, VirtualPosition};
use noodles::core::{Position, Region};
use noodles::csi::{self as csi, BinningIndex};
use noodles::tabix;

#[derive(Debug)]
pub struct QueryResult(pub Vec<u32>, pub Vec<u32>, pub u32);

pub enum CallableSites {
    D4(D4CallableSites),
    Bed(BedCallableSites),
}

impl CallableSites {
    pub fn query(
        &mut self,
        chrom: &str,
        begin: u32,
        end: u32,
        ploidy: u32,
        skip_sites: &HashSet<u32>,
    ) -> Result<QueryResult> {
        match self {
            CallableSites::D4(reader) => reader.query(chrom, begin - 1, end, ploidy, skip_sites),
            CallableSites::Bed(reader) => reader.query(chrom, begin, end, ploidy, skip_sites),
        }
    }
}

pub struct BedCallableSites {
    decoder: bgzf::Reader<File>,
    index: csi::binning_index::Index<Vec<VirtualPosition>>,
    contigs: IndexSet<String>,
    num_samples: Vec<usize>,
}

impl BedCallableSites {
    pub fn from_file<P: AsRef<Path>>(src: P, num_samples: Vec<usize>) -> Result<Self> {
        let index_src = format!("{}.tbi", src.as_ref().display());
        let index =
            tabix::read(&index_src).context(format!("Couldn't open index: {:?}", &index_src))?;

        let header = index.header().expect("missing tabix header");
        let contigs = header
            .reference_sequence_names()
            .into_iter()
            .map(|bs| String::from_utf8_lossy(&bs).into_owned())
            .collect();

        let decoder = bgzf::reader::Builder.build_from_path(src)?;

        Ok(BedCallableSites {
            decoder,
            index,
            contigs,
            num_samples,
        })
    }

    pub fn query(
        &mut self,
        chrom: &str,
        begin: u32,
        end: u32,
        ploidy: u32,
        skip_sites: &HashSet<u32>,
    ) -> Result<QueryResult> {
        let begin_pos = Position::new(begin as usize).context("Invalid start position")?;
        let end_pos = Position::new(end as usize).context("Invalid end position")?;
        let contig_id = self
            .contigs
            .get_index_of(chrom)
            .context("Invalid contig name")?;
        let region = Region::new(BString::from(chrom), begin_pos..=end_pos);
        let chunks = self.index.query(contig_id, region.interval())?;
        let query = csi::io::Query::new(&mut self.decoder, chunks);

        let mut reader = bed::io::Reader::<3, _>::new(query);
        let mut record = bed::Record::default();

        let mut within_comps = vec![0; self.num_samples.len()];
        let mut dxy_comps = if self.num_samples.len() > 1 {
            vec![0; self.num_samples.len() * (self.num_samples.len() - 1) / 2]
        } else {
            vec![0]
        };
        let mut callable_sites = 0;
        while reader.read_record(&mut record)? != 0 {
            let record_start = record.feature_start()?;
            let record_end = record.feature_end().expect("missing feature end")?;
            let record_interval = (record_start..=record_end).into();

            if !region.interval().intersects(record_interval) {
                continue;
            }

            let overlap_start = std::cmp::max(record_start.get(), begin as usize);
            let overlap_end = std::cmp::min(record_end.get(), end as usize);
            let sites = if skip_sites.is_empty() {
                overlap_end - overlap_start + 1
            } else {
                (overlap_start..=overlap_end)
                    .filter(|&pos| !skip_sites.contains(&(pos as u32)))
                    .count()
            };
            callable_sites += sites as u32;
            for (pop_idx, &num_samples) in self.num_samples.iter().enumerate() {
                let haps = ploidy * num_samples as u32;
                trace!(
                    "Popidx: {} haps: {}, sites: {}, olstart: {}, olend: {}",
                    pop_idx,
                    haps,
                    sites,
                    overlap_start,
                    overlap_end
                );
                let within_comp = if haps > 1 {
                    (haps as u64 * (haps as u64 - 1) / 2) as u32
                } else {
                    0
                };
                within_comps[pop_idx] += within_comp * sites as u32;
            }

            if self.num_samples.len() > 1 {
                let mut comb_idx = 0;
                for i in 0..self.num_samples.len() {
                    let haps_i = ploidy * self.num_samples[i] as u32;
                    for j in (i + 1)..self.num_samples.len() {
                        let haps_j = ploidy * self.num_samples[j] as u32;
                        dxy_comps[comb_idx] += (haps_i * haps_j) * sites as u32;
                        comb_idx += 1;
                    }
                }
            }
        }

        Ok(QueryResult(within_comps, dxy_comps, callable_sites))
    }
}
pub struct D4CallableSites {
    readers: Vec<D4TrackReader<File>>,
}

impl D4CallableSites {
    pub fn from_file<P: AsRef<Path>>(pops: Vec<&str>, d4_file: P) -> Result<Self> {
        let d4_file_path = d4_file.as_ref();
        let mut readers: Vec<D4TrackReader<File>> = Vec::new();

        if pops.len() > 1 {
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

    pub fn query(
        &mut self,
        chrom: &str,
        begin: u32,
        end: u32,
        ploidy: u32,
        skip_sites: &HashSet<u32>,
    ) -> Result<QueryResult> {
        let num_pops = self.readers.len();
        trace!("querying d4 {}:{}-{}", chrom, begin, end);
        let mut views: Vec<D4TrackView<File>> = self
            .readers
            .iter_mut()
            .map(|x| x.get_view(chrom, begin, end).unwrap())
            .collect();

        let mut within_comps: Vec<u32> = vec![0; num_pops];
        let mut callable_sites = 0;
        let num_pop_combs = if num_pops > 1 {
            num_pops * (num_pops - 1) / 2
        } else {
            0
        };

        let mut dxy_comps = if num_pop_combs > 0 {
            vec![0; num_pop_combs]
        } else {
            vec![0]
        };

        for pos in begin..end {
            // Get the callable individuals (values) for each population at this position

            let mut any_callable = false;
            let values: Vec<u32> = views
                .iter_mut()
                .map(|view| {
                    let (reported_pos, value) = view.next().unwrap().unwrap();

                    assert_eq!(reported_pos, pos);
                    let value = value as u32;
                    if value != 0 {
                        any_callable = true;
                    }

                    value
                })
                .collect();

            if skip_sites.contains(&(pos + 1)) || !any_callable {
                continue;
            }

            callable_sites += 1;

            // Calculate within-population comparisons
            for (pop_idx, &callable_indvs) in values.iter().enumerate() {
                let haps = ploidy * callable_indvs;

                let within_comp = if haps > 1 {
                    (haps as u64 * (haps as u64 - 1) / 2) as u32
                } else {
                    0
                };
                within_comps[pop_idx] += within_comp;
            }

            // Calculate between-population (dxy) comparisons
            if num_pop_combs > 0 {
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
            "Callable Sites - {}:{}-{} {:?}, {}",
            chrom, begin, end, within_comps, callable_sites
        );
        Ok(QueryResult(within_comps, dxy_comps, callable_sites))
    }
}

#[cfg(test)]
mod tests {
    use std::collections::HashSet;
    use std::fs::File;
    use std::path::{Path, PathBuf};

    use serde::Deserialize;

    use super::*;

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
        let mut d4callable = D4CallableSites::from_file(Vec::<&str>::default(), fp)?;

        let QueryResult(within_comps, _, _) = d4callable.query(
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
        let mut d4callable = D4CallableSites::from_file(Vec::<&str>::default(), fp)?;

        for result in reader.deserialize() {
            let record: NoPopsTruthRecord = result?;
            let QueryResult(within_comps, _, _) = d4callable.query(
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
        let mut d4callable = D4CallableSites::from_file(Vec::<&str>::default(), fp)?;
        let mut skip_sites = HashSet::new();
        skip_sites.insert(10 as u32);
        for result in reader.deserialize() {
            let mut record: NoPopsTruthRecord = result?;
            let QueryResult(within_comps, _, _) =
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
        let mut d4callable = D4CallableSites::from_file(pops, fp)?;

        for result in reader.deserialize() {
            let record: PopsTruthRecord = result?;
            let QueryResult(within_comps, dxy_comps, _) = d4callable.query(
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

            assert_eq!(
                dxy_comps[0], record.dxy,
                "Mismatch in dxy for range {}:{}-{}",
                record.chrom, record.start, record.end
            );
        }
        Ok(())
    }

    #[test]
    fn test_bed_query() -> Result<()> {
        init();
        let test_path = PathBuf::from("tests/data/stat/callable_sites.bed.gz");
        let num_samples = vec![5, 3]; // Two populations with 5 and 3 samples
        let mut reader = BedCallableSites::from_file(test_path, num_samples)?;

        // Test chr1:50-150 (overlaps first region)
        let QueryResult(within, dxy, _) = reader.query("chr1", 50, 150, 2, &HashSet::new())?;
        assert_eq!(within, vec![2295, 765]); // pop1 10 haplotypes = 10c2 = 45 * 100 =

        Ok(())
    }
}
