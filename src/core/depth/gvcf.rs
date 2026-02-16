use color_eyre::{eyre::WrapErr, Result};
use log::{debug, trace, warn};
use noodles::bgzf;
use noodles::core::Region;
use noodles::tabix;
use noodles::vcf;
use noodles::vcf::variant::record::info::field::Value::Integer as info_int;
use noodles::vcf::variant::record::samples::series::Value::Integer as samples_int;
use noodles::vcf::variant::record::samples::Series;
use noodles::vcf::Header;
use std::fs::File;
use std::path::{Path, PathBuf};

use super::DepthSource;
use crate::core::contig::{Contig, ContigSet};

pub struct GvcfReader {
    inner: vcf::io::Reader<bgzf::Reader<File>>,
    index: tabix::Index,
    header: Header,
    src: PathBuf,
    sample_name: String,
    min_gq: isize,
}

impl GvcfReader {
    pub fn new<P: AsRef<Path>>(src: P, sample_name: &str, min_gq: Option<isize>) -> Result<Self> {
        let src_path = src.as_ref().to_path_buf();

        let index = tabix::fs::read(format!("{}.tbi", src_path.display()))
            .wrap_err_with(|| format!("Failed to read tabix index for: {}", src_path.display()))?;

        let file = File::open(&src_path)
            .wrap_err_with(|| format!("Failed to open GVCF file: {}", src_path.display()))?;

        let mut reader = vcf::io::Reader::new(bgzf::Reader::new(file));

        let header = reader
            .read_header()
            .wrap_err_with(|| format!("Failed to read GVCF header: {}", src_path.display()))?;

        Ok(Self {
            inner: reader,
            index,
            header,
            src: src_path,
            sample_name: sample_name.to_string(),
            min_gq: min_gq.unwrap_or(1),
        })
    }
}

impl DepthSource for GvcfReader {
    fn read_depths(&mut self, chrom: &str, start: u32, end: u32) -> Result<Vec<u32>> {
        let len = (end - start) as usize;
        let mut depths = vec![0u32; len];

        let region: Region = format!("{}:{}-{}", chrom, start + 1, end)
            .parse()
            .wrap_err_with(|| format!("Invalid region: {}:{}-{}", chrom, start, end))?;

        let query = self
            .inner
            .query(&self.header, &self.index, &region)
            .wrap_err_with(|| format!("Failed to query region: {}:{}-{}", chrom, start, end))?;

        trace!(
            "Querying GVCF {} region {}:{}-{}",
            self.src.display(),
            chrom,
            start,
            end
        );

        for result in query {
            let record = result.wrap_err("Failed to read GVCF record")?;

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
            let Some(format_dp) = samples.select("DP") else {
                debug!(
                    "Missing DP format field in file {}\nOffending record: {:?}",
                    self.src.display(),
                    record
                );
                continue;
            };

            let Some(Some(Ok(samples_int(dp)))) = format_dp.get(&self.header, 0) else {
                debug!(
                    "Missing DP value in file {}\nOffending record: {:?}",
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
                    "Missing GQ value in file {}\nOffending record: {:?}",
                    self.src.display(),
                    record
                );
                continue;
            };

            let array_start = startpos.get() - 1;
            let array_end = endpos as usize;

            trace!(
                "Record at {}:{}-{} DP={}, GQ={}",
                chrom,
                array_start,
                array_end,
                dp,
                gq
            );

            if gq >= self.min_gq as i32 {
                for pos in array_start..array_end {
                    let idx = pos.saturating_sub(start as usize);
                    if idx < depths.len() {
                        depths[idx] = dp as u32;
                    }
                }
            }
        }

        Ok(depths)
    }

    fn contigs(&self) -> ContigSet {
        let contigs: Vec<Contig> = self
            .header
            .contigs()
            .iter()
            .filter_map(|(name, contig_map)| {
                contig_map
                    .length()
                    .map(|len| Contig::new(name.to_string(), len))
            })
            .collect();

        ContigSet::new(contigs)
    }

    fn sample_name(&self) -> &str {
        &self.sample_name
    }
}