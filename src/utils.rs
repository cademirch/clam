use anyhow::{bail, Result};
use fnv::FnvHashMap;
use indexmap::IndexSet;
use noodles::vcf::{self};
use noodles::{bed, core::Region};
use serde::Deserialize;
use std::collections::HashSet;
use std::fs::File;
use std::path::Path;
use camino::Utf8PathBuf;
use std::io::{BufReader, BufRead};


#[derive(Debug, Deserialize)]
struct PopFileRecord {
    sample: String,
    population_name: String,
}

pub fn count_combinations(n: u32, r: u32) -> u32 {
    if r > n {
        0
    } else {
        (1..=r).fold(1, |acc, val| acc * (n - val + 1) / val)
    }
}

#[derive(Debug, Clone)]
pub struct PopulationMapping {
    pub pop_idx_to_pop_name: Option<Vec<String>>,
    pub pop_idx_to_sample_names: Vec<Vec<String>>,
    pub num_populations: usize,
    pub samp_idx_to_pop_idx: FnvHashMap<usize, usize>,
    pub sample_name_to_sample_idx: FnvHashMap<String, usize>
}

impl PopulationMapping {
    pub fn default(header: &vcf::Header) -> Self {
        let mut sample_to_pop_idx = FnvHashMap::default();
        let mut pop_idx_to_sample_names = vec![Vec::new()];
        let mut sample_name_to_sample_idx = FnvHashMap::default(); // Initialize the new map
        let pop_idx_to_pop_name = None;

        for (idx, sample_name) in header.sample_names().iter().enumerate() {
            let sample_name_str = sample_name.clone(); // Clone sample name
            sample_to_pop_idx.insert(idx, 0); // Map all samples to population 0
            pop_idx_to_sample_names[0].push(sample_name_str.clone()); // Add all samples to Population_0
            sample_name_to_sample_idx.insert(sample_name_str, idx); // Map sample name to its index
        }

        PopulationMapping {
            pop_idx_to_pop_name,
            pop_idx_to_sample_names,
            num_populations: 1,
            samp_idx_to_pop_idx: sample_to_pop_idx,
            sample_name_to_sample_idx, // Add the new map
        }
    }


    /// Creates a PopulationMapping from a file containing sample and population data.
    pub fn from_path<P: AsRef<Path>>(
        pop_file: P,
        header_samples: Option<&IndexSet<String>>, // Make header_samples optional
    ) -> Result<Self> {
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
            let pop_idx = if let Some(idx) = pop_idx_to_name.iter().position(|p| p == population_name) {
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
    
        // Optional creation of mappings dependent on header_samples
        let mut samp_idx_to_pop_idx = FnvHashMap::default();
        let mut sample_name_to_sample_idx = FnvHashMap::default();
        if let Some(header_samples) = header_samples {
            for (pop_idx, sample_names) in pop_idx_to_sample_names.iter().enumerate() {
                for sample in sample_names.iter() {
                    if let Some(sample_idx) = header_samples.get_index_of(sample) {
                        samp_idx_to_pop_idx.insert(sample_idx, pop_idx);
                        sample_name_to_sample_idx.insert(sample.clone(), sample_idx);
                    } else {
                        bail!("Sample: '{}' not found in VCF header!", sample);
                    }
                }
            }
        }
    
        Ok(PopulationMapping {
            pop_idx_to_pop_name: Some(pop_idx_to_name),
            pop_idx_to_sample_names,
            num_populations,
            samp_idx_to_pop_idx,
            sample_name_to_sample_idx,
        })
    }

    pub fn get_sample_refs(&self, pop_idx: usize) -> Vec<&str> {
        self.pop_idx_to_sample_names[pop_idx]
            .iter()
            .map(|s| s.as_str())
            .collect()
    }

    pub fn get_popname_refs(&self) -> Option<Vec<&str>> {
        self.pop_idx_to_pop_name
            .as_ref()
            .map(|vec| vec.iter().map(String::as_str).collect())
    }
}


pub fn get_exclude_chromosomes(
    exclude: &Option<Vec<String>>,
    exclude_file: &Option<Utf8PathBuf>,
) -> Result<HashSet<String>> {
    if let Some(file_path) = exclude_file {
        // Read chromosomes from the file and collect them into a HashSet
        let file = File::open(file_path)?;
        let reader = BufReader::new(file);
        let exclude_set: HashSet<String> = reader
            .lines()
            .filter_map(Result::ok) // Ignore any lines that fail to read
            .collect();
        Ok(exclude_set)
    } else if let Some(chromosomes) = exclude {
        // Convert Vec<String> to HashSet<String> directly
        Ok(chromosomes.iter().cloned().collect())
    } else {
        // If neither option is provided, return an empty HashSet
        Ok(HashSet::new())
    }
}

pub fn read_bed_regions<P: AsRef<Path>>(bed_file: P) -> Result<Vec<Region>> {
    let mut res = vec![];
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
        let start = match record.feature_start() {
            Ok(pos) => pos,
            Err(e) => {
                bail!("Bad bed record start: {:?}. Error: {:?}", record, e);
            }
        };
        let end = match record.feature_end() {
            Some(Ok(pos)) => pos,
            Some(Err(e)) => bail!("Bad bed record end: {:?}. Error: {:?}", record, e),
            None => bail!("Bed record missing end: {:?}", record),
        };
        let region = Region::new(chrom, start..=end);
        
        res.push(region);
    }


    Ok(res)
}

