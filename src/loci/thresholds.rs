use std::collections::hash_map::Entry;
use std::collections::HashMap;
use std::fs::File;
use std::path::Path;

use anyhow::Result;
use log::warn;
use rayon::prelude::*;
use serde::Deserialize;

#[derive(Clone)]
pub enum Thresholds {
    Fixed((f64, f64)),                          // For single threshold tuple
    PerChromosome(HashMap<String, (f64, f64)>), // For chromosome-specific thresholds
}

#[derive(Debug, Deserialize)]
struct ThresholdRecord {
    chrom: String,
    min: f64,
    max: f64,
}

pub fn read_threshold_file<P: AsRef<Path>>(
    threshold_file: P,
) -> Result<HashMap<String, (f64, f64)>> {
    let file = File::open(&threshold_file).expect(&format!(
        "Failed to open threshold file: {}",
        threshold_file.as_ref().display()
    ));
    let mut reader = csv::ReaderBuilder::new().delimiter(b'\t').from_reader(file);
    let mut filter_map: HashMap<String, (f64, f64)> = HashMap::new();

    for result in reader.deserialize() {
        let record: ThresholdRecord = result?;

        match filter_map.entry(record.chrom) {
            Entry::Vacant(entry) => {
                entry.insert((record.min, record.max));
            }
            Entry::Occupied(mut entry) => {
                warn!(
                    "Duplicate chromosome in thresholds file! {} already exists with min: {}, max: {}. Overwriting with new values min: {}, max: {}.",
                    entry.key(),
                    entry.get().0,
                    entry.get().1,
                    record.min,
                    record.max
                );
                entry.insert((record.min, record.max));
            }
        }
    }

    Ok(filter_map)
}
