use color_eyre::Result;
use indicatif::{ProgressBar, ProgressStyle};
use noodles::csi::{self, BinningIndex};
use noodles::tabix;
use serde::Deserialize;
use std::collections::HashMap;
use std::ffi::OsString;
use std::io;
use std::path::{Path, PathBuf};

/// Read the index of a bgzipped, coordinate-sorted file: `<src>.tbi` if
/// present, otherwise `<src>.csi`. CSI is needed for references with
/// sequences longer than 2^29 bp, which tabix cannot index. Mirrors how
/// noodles' VCF indexed reader finds its index.
pub fn read_tabix_or_csi_index(src: &Path) -> io::Result<Box<dyn BinningIndex + Send + Sync>> {
    match tabix::fs::read(index_path(src, "tbi")) {
        Ok(index) => Ok(Box::new(index)),
        Err(e) if e.kind() == io::ErrorKind::NotFound => {
            let index = csi::fs::read(index_path(src, "csi"))?;
            Ok(Box::new(index))
        }
        Err(e) => Err(e),
    }
}

fn index_path(src: &Path, ext: &str) -> PathBuf {
    let mut s = OsString::from(src);
    s.push(".");
    s.push(ext);
    PathBuf::from(s)
}

// Helper to create a consistent spinner
pub fn create_spinner(message: &str) -> ProgressBar {
    let pb = ProgressBar::new_spinner();
    pb.set_style(
        ProgressStyle::default_spinner()
            .template("{spinner:.green} {msg}")
            .unwrap(),
    );
    pb.set_message(message.to_string());
    pb.enable_steady_tick(std::time::Duration::from_millis(100));
    pb
}

/// Helper to create a consistent progress bar
pub fn create_progress_bar(total: usize) -> ProgressBar {
    let pb = ProgressBar::new(total as u64);
    pb.set_style(
        ProgressStyle::default_bar()
            .template("{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {pos}/{len} chunks ({per_sec}, {eta})")
            .unwrap()
            .progress_chars("#>-")
    );
    pb
}

#[derive(Debug, Deserialize)]
struct ThresholdRecord {
    contig: String,
    min_depth: f64,
    max_depth: f64,
}

pub fn parse_contig_thresholds(path: &Path) -> Result<HashMap<String, (f64, f64)>> {
    let mut rdr = csv::ReaderBuilder::new().delimiter(b'\t').from_path(path)?;

    let mut thresholds = HashMap::new();
    for result in rdr.deserialize() {
        let record: ThresholdRecord = result?;
        thresholds.insert(record.contig, (record.min_depth, record.max_depth));
    }

    Ok(thresholds)
}
