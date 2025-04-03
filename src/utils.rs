use std::collections::hash_map::Entry;
use std::collections::HashSet;
use std::fs::File;
use std::io::{BufRead, BufReader, Read};
use std::path::Path;

use anyhow::{bail, Context, Result};
use camino::Utf8PathBuf;
use fnv::FnvHashMap;
use indexmap::{indexset, IndexSet};
use log::warn;
use noodles::bed;
use noodles::core::Region;
use serde::Deserialize;

const DEFAULT_POPULATION_NAME: &str = "default_population";

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
    population_names: IndexSet<String>, // index corresponds to population_idx
    sample_names: Vec<Vec<String>>, 
    sample_name_lookup: FnvHashMap<String, (usize, usize)>, // key is sample name. value is (population_idx, sample's idx within population (internal_idx))
    sample_idx_lookup: Option<FnvHashMap<usize, (usize, usize)>>, // key is sample idx in vcf (global_idx). value is (population_idx, sample's idx within population (internal_idx))
}

impl PopulationMapping {
    /// No population file, just sample names from VCF header. Makes all sample as one population.
    pub fn default(vcf_sample_names: &IndexSet<String>) -> Self {
        let population_names = indexset! {DEFAULT_POPULATION_NAME.to_string()};
        let mut sample_name_lookup = FnvHashMap::default();
        let mut sample_idx_lookup = FnvHashMap::default();
        let mut sample_names: Vec<Vec<String>> = vec![Vec::default();1];
        for (global_idx, sample_name) in vcf_sample_names.iter().enumerate() {
            let population_idx = 0; // Default population
            let internal_idx = global_idx;
            sample_names[0].push(sample_name.clone());
            sample_name_lookup.insert(sample_name.clone(), (population_idx, internal_idx));
            sample_idx_lookup.insert(global_idx, (population_idx, internal_idx));
        }

        Self {
            population_names,
            sample_names,
            sample_name_lookup,
            sample_idx_lookup: Some(sample_idx_lookup),
        }
    }

    /// Create from population file and optional vcf header (for loci support)
    pub fn from_path<R: Read>(
        reader: R,
        vcf_sample_names: Option<&IndexSet<String>>, // from vcf header
    ) -> Result<Self> {
        let mut csv_reader = csv::ReaderBuilder::new()
            .delimiter(b'\t')
            .has_headers(false)
            .from_reader(reader);

        let mut population_names = IndexSet::default();
        let mut sample_name_lookup = FnvHashMap::default();
        let mut sample_idx_lookup = if vcf_sample_names.is_some() {
            Some(FnvHashMap::default())
        } else {
            None
        };

        let mut sample_names: Vec<Vec<String>> = Vec::default();

        for result in csv_reader.deserialize() {
            let record: PopFileRecord = result?;

            let new_population = population_names.insert(record.population_name.clone());
            let population_idx = population_names
                .get_index_of(&record.population_name)
                .unwrap();

            if new_population {
                sample_names.push(vec![record.sample.to_string()])
            } else {
                sample_names[population_idx].push(record.sample.to_string())
            };
        }

        for (population_idx, samples) in sample_names.iter().enumerate() {
            for (internal_idx, sample_name) in samples.iter().enumerate() {
                match sample_name_lookup.entry(sample_name.clone()) {
                    Entry::Vacant(entry) => {
                        entry.insert((population_idx, internal_idx));
                    }
                    Entry::Occupied(entry) => {
                        bail!("Duplicated sample name: {} in population file!", {
                            entry.key()
                        })
                    }
                }

                if let Some(ref mut sample_idx_lookup) = sample_idx_lookup {
                    let global_idx = vcf_sample_names
                        .unwrap()
                        .get_index_of(&sample_name.to_string())
                        .context(format!(
                            "Sample: '{}' not in VCF header!",
                            &sample_name.to_string()
                        ))?;
                    sample_idx_lookup.insert(global_idx, (population_idx, internal_idx));
                }
            }
        }

        Ok(Self {
            population_names,
            sample_names,
            sample_name_lookup,
            sample_idx_lookup,
        })
    }
    
    pub fn get_samples_per_population(&self) -> Vec<Vec<&str>> {
        self.sample_names
            .iter()
            .map(|population| population.iter().map(String::as_str).collect())
            .collect()
    }

    pub fn get_sample_counts_per_population(&self) -> Vec<usize> {
        let counts: Vec<usize> = self.sample_names.iter().map(|x| x.len()).collect();
        counts
    }
    
    pub fn num_populations(&self) -> usize {
        self.population_names.len()
    }
    pub fn lookup_sample_name(&self, sample_name: &str) -> Result<(usize, usize)> {
        let res = self
            .sample_name_lookup
            .get(sample_name)
            .context("Sample name doesn't exist!")?;
        Ok(*res)
    }

    pub fn lookup_sample_idx(&self, sample_idx: usize) -> Result<(usize, usize)> {
        let res = if let Some(ref sample_idx_lookup) = self.sample_idx_lookup {
            Ok(*sample_idx_lookup
                .get(&sample_idx)
                .context("Sample idx doesn't exist!")?)
        } else {
            bail!("This shouldn't happen! Please open a github issue.")
        };

        res
    }

    pub fn get_sample_refs(&self, pop_idx: usize) -> Vec<&str> {
        self.sample_name_lookup
            .iter()
            .filter_map(|(sample_name, &(population_idx, _))| {
                if population_idx == pop_idx {
                    Some(sample_name.as_str())
                } else {
                    None
                }
            })
            .collect()
    }

    pub fn get_popname_refs(&self) -> Vec<&str> {
        self.population_names.iter().map(String::as_str).collect()
    }
    pub fn validate_sample_coverage(&self, vcf_sample_names: &IndexSet<String>) -> Result<()> {
        // Collect all sample names from the population file
        let pop_file_samples: IndexSet<_> = self.sample_name_lookup.keys().cloned().collect();

        // Find samples in the population file but not in the VCF header
        let missing_in_vcf: IndexSet<_> = pop_file_samples.difference(vcf_sample_names).collect();

        // Find samples in the VCF header but not in the population file
        let missing_in_pop_file: IndexSet<_> =
            vcf_sample_names.difference(&pop_file_samples).collect();

        if !missing_in_vcf.is_empty() {
            bail!(
                "The following samples are in the population file but missing in the VCF header: {:?}",
                missing_in_vcf
            );
        }

        if !missing_in_pop_file.is_empty() {
            warn!(
                "The following samples are in the VCF header but missing in the population file: {:?}\nThey will be skipped!",
                missing_in_pop_file
            );
        }

        
        Ok(())
    }
}

pub fn get_exclude_chromosomes(
    exclude: &Option<Vec<String>>,
    exclude_file: &Option<Utf8PathBuf>,
) -> Result<Option<HashSet<String>>> {
    if let Some(file_path) = exclude_file {
        // Read chromosomes from the file and collect them into a HashSet
        let file = File::open(file_path)?;
        let reader = BufReader::new(file);
        let exclude_set: HashSet<String> = reader
            .lines()
            .filter_map(Result::ok) // Ignore any lines that fail to read
            .collect();
        Ok(Some(exclude_set))
    } else if let Some(chromosomes) = exclude {
        // Convert Vec<String> to HashSet<String> directly
        Ok(Some(chromosomes.iter().cloned().collect()))
    } else {
        // If neither option is provided, return an empty HashSet
        Ok(None)
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
        Ok(_) => true,  // continue loop
        Err(e) => {
            bail!(
                "Error reading record: {}. Offending record: {:?}",
                e,
                record
            );
        }
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
