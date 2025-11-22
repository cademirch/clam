use color_eyre::{eyre::WrapErr, Result};
use log::{debug, trace, warn};
use noodles::bgzf::Reader;
use noodles::core::Region;
use noodles::vcf::io::IndexedReader;
use noodles::vcf::variant::record::info::field::Value::Integer as info_int;
use noodles::vcf::variant::record::samples::series::Value::Integer as samples_int;
use noodles::vcf::variant::record::samples::Series;
use noodles::vcf::Header;
use std::fs::File;
use std::path::{Path, PathBuf};

use super::DepthSource;
use crate::core::contig::{Contig, ContigSet};

pub struct GvcfReader {
    inner: IndexedReader<Reader<File>>,
    header: Header,
    src: PathBuf,
    sample_name: String,
}

impl GvcfReader {
    pub fn new<P: AsRef<Path>>(src: P, sample_name: &str) -> Result<Self> {
        let src_path = src.as_ref().to_path_buf();

        let mut reader = noodles::vcf::io::indexed_reader::Builder::default()
            .build_from_path(&src_path)
            .wrap_err_with(|| format!("Failed to open GVCF file: {}", src_path.display()))?;

        let header = reader
            .read_header()
            .wrap_err_with(|| format!("Failed to read GVCF file header: {}", src_path.display()))?;

        Ok(Self {
            inner: reader,
            header,
            src: src_path,
            sample_name: sample_name.to_string(),
        })
    }
}

impl DepthSource for GvcfReader {
    fn read_depths(&mut self, chrom: &str, start: u32, end: u32) -> Result<Vec<u32>> {
        let len = (end - start) as usize;
        let mut depths = vec![0u32; len];

        // Create region query (noodles uses 1-based coordinates)
        let region: Region = format!("{}:{}-{}", chrom, start + 1, end)
            .parse()
            .wrap_err_with(|| format!("Invalid region: {}:{}-{}", chrom, start, end))?;

        let query = self
            .inner
            .query(&self.header, &region)
            .wrap_err_with(|| format!("Failed to query region: {}:{}-{}", chrom, start, end))?;

        trace!("Querying GVCF region {}:{}-{}", chrom, start, end);

        for result in query {
            let record = result.wrap_err("Failed to read GVCF record")?;

            // Get start position
            let Some(Ok(startpos)) = record.variant_start() else {
                warn!(
                    "Invalid start position in file {}\nOffending record: {:?}",
                    self.src.display(),
                    record
                );
                continue;
            };

            // Get END position from INFO field (or default to start + 1)
            let info = record.info();
            let endpos = match info.get(&self.header, "END") {
                Some(Ok(Some(info_int(end)))) => end,
                _ => (startpos.get() + 1) as i32,
            };

            // Get DP (depth) value
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

            // Get GQ (genotype quality) value
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

            // Convert to 0-based array indices
            let array_start = (startpos.get() - 1) as usize;
            let array_end = endpos as usize;

            trace!(
                "Record at {}:{}-{} DP={}, GQ={}",
                chrom,
                array_start,
                array_end,
                dp,
                gq
            );

            // Fill in depths for this interval
            // If GQ is 0, treat depth as 0 (not callable)
            let depth_value = if gq > 0 { dp as u32 } else { 0 };

            for pos in array_start..array_end {
                let idx = pos.saturating_sub(start as usize);
                if idx < depths.len() {
                    depths[idx] = depth_value;
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
                    .map(|len| Contig::new(name.to_string(), len as usize))
            })
            .collect();

        ContigSet::new(contigs)
    }

    fn sample_name(&self) -> &str {
        &self.sample_name
    }
}

