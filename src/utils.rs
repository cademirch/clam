use anyhow::{bail, Result};
use fnv::FnvHashMap;
use log::warn;
use noodles::{bed, core::Region};
use serde::Deserialize;
use std::collections::{hash_map::Entry, HashMap};
use std::fs::File;
use std::path::Path;

#[derive(Debug, Clone)]
pub struct PopulationMapping {
    pub sample_name_to_pop_idx: FnvHashMap<String, usize>,
    pub pop_idx_to_pop_name: Vec<String>,
    pub pop_idx_to_sample_names: Vec<Vec<String>>,
    pub num_populations: usize,
}
#[derive(Debug, Deserialize)]
struct PopFileRecord {
    sample: String,
    population_name: String,
}

impl PopulationMapping {
    /// Creates a PopulationMapping from a file containing sample and population data.
    pub fn from_path<P: AsRef<Path>>(pop_file: P) -> Result<PopulationMapping> {
        let file = File::open(&pop_file).expect(&format!(
            "Failed to open population file: {}",
            pop_file.as_ref().display()
        ));
        let mut reader = csv::ReaderBuilder::new().delimiter(b'\t').from_reader(file);

        let mut sample_name_to_pop_idx = FnvHashMap::default();
        let mut pop_idx_to_name = Vec::new();
        let mut sample_records: Vec<(String, String)> = Vec::new();

        for result in reader.deserialize() {
            let record: PopFileRecord = result?;
            sample_records.push((record.sample, record.population_name));
        }

        for (sample_name, population_name) in &sample_records {
            let pop_idx =
                if let Some(idx) = pop_idx_to_name.iter().position(|p| p == population_name) {
                    idx
                } else {
                    let idx = pop_idx_to_name.len();
                    pop_idx_to_name.push(population_name.clone());
                    idx
                };

            sample_name_to_pop_idx.insert(sample_name.clone(), pop_idx);
        }

        let num_populations = pop_idx_to_name.len();

        let mut pop_idx_to_sample_names: Vec<Vec<String>> = vec![Vec::new(); num_populations];
        for (sample_name, pop_idx) in sample_name_to_pop_idx.iter() {
            pop_idx_to_sample_names[*pop_idx].push(sample_name.clone());
        }

        Ok(PopulationMapping {
            sample_name_to_pop_idx,
            pop_idx_to_pop_name: pop_idx_to_name,
            pop_idx_to_sample_names,
            num_populations,
        })
    }

    pub fn get_sample_refs(&self, pop_idx: usize) -> Vec<&str> {
        self.pop_idx_to_sample_names[pop_idx]
            .iter()
            .map(|s| s.as_str())
            .collect()
    }
    pub fn get_popname_refs(&self) -> Vec<&str> {
        self.pop_idx_to_pop_name
            .iter()
            .map(|s| s.as_str())
            .collect()
    }
}

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
        "Failed to open population file: {}",
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

pub fn read_bed_regions<P: AsRef<Path>>(bed_file: P) -> Result<HashMap<String, Vec<Region>>> {
    let mut res: HashMap<String, Vec<Region>> = HashMap::new();
    let mut rdr = bed::io::reader::Builder::<3>::default()
        .build_from_path(bed_file.as_ref())
        .expect(&format!(
            "Failed to open bed file: {}",
            bed_file.as_ref().display()
        ));
    let mut record = bed::Record::<3>::default();
    
    while match rdr.read_record(&mut record) {
        Ok(0) => false, // no more records break loop
        Ok(_) => true, // continue loop
        Err(e) => {
            bail!("Error reading record: {}. Offending record: {:?}", e, record);
            
        },
    } {
        let chrom = record.reference_sequence_name();
        let start = record.feature_start().expect(&format!("Bad bed record start: {:?}", record));
        let end = match record.feature_end() {
            Some(pos) => pos,
            None => bail!("Bed record missing end: {:?}", record)
        }?;
        let region = Region::new(chrom, start..=end);
        
        match res.entry(chrom.to_string()) {
            Entry::Vacant(e) => {e.insert(vec![region]);},
            Entry::Occupied(mut e) => {e.get_mut().push(region);}
        };
    }


    Ok(res)
}
