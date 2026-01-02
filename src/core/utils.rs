use color_eyre::Result;
use indicatif::{ProgressBar, ProgressStyle};
use serde::Deserialize;
use std::collections::HashMap;
use std::path::Path;

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
