use color_eyre::{
    eyre::{bail, Context, ContextCompat},
    Result,
};
use indexmap::IndexMap;
use std::collections::{HashMap, HashSet};
use std::fmt;
use std::fs::File;
use std::path::Path;
/// A single population with its samples
#[derive(Debug, Clone)]
pub struct Population {
    pub name: String,
    samples: Vec<String>,
}

impl Population {
    pub fn new(name: String, samples: Vec<String>) -> Self {
        Self { name, samples }
    }

    pub fn samples(&self) -> &[String] {
        &self.samples
    }

    pub fn num_samples(&self) -> usize {
        self.samples.len()
    }
}

impl fmt::Display for Population {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.name)
    }
}
/// Maps samples to populations
#[derive(Debug, Clone)]
pub struct PopulationMap {
    populations: IndexMap<String, Population>,
    // sample_name -> (population_index, sample_index_within_population)
    sample_lookup: HashMap<String, (usize, usize)>,
}

impl PopulationMap {
    /// Create from a tab-separated file: sample\tpopulation
    pub fn from_file(path: impl AsRef<Path>) -> Result<Self> {
        let mut populations: IndexMap<String, Vec<String>> = IndexMap::new();

        let file = File::open(path.as_ref()).wrap_err_with(|| {
            format!(
                "Failed to open population file: {}",
                path.as_ref().display()
            )
        })?;

        let mut reader = csv::ReaderBuilder::new()
            .delimiter(b'\t')
            .has_headers(false)
            .from_reader(file);

        for (line_num, result) in reader.records().enumerate() {
            let record = result.wrap_err_with(|| {
                format!("Failed to parse line {} in population file", line_num + 1)
            })?;

            let sample_name = record
                .get(0)
                .wrap_err_with(|| format!("Missing sample name on line {}", line_num + 1))?
                .to_string();

            let pop_name = record
                .get(1)
                .wrap_err_with(|| format!("Missing population name on line {}", line_num + 1))?
                .to_string();

            populations.entry(pop_name).or_default().push(sample_name);
        }

        Self::from_populations(populations)
    }

    /// Create default single population from a list of samples
    pub fn default_from_samples(samples: Vec<String>) -> Self {
        let mut populations = IndexMap::new();
        populations.insert("default".to_string(), samples);
        Self::from_populations(populations).expect("Default population should always be valid")
    }

    /// Internal constructor that builds lookup tables
    fn from_populations(pop_data: IndexMap<String, Vec<String>>) -> Result<Self> {
        let mut populations = IndexMap::new();
        let mut sample_lookup = HashMap::new();

        for (pop_idx, (pop_name, samples)) in pop_data.into_iter().enumerate() {
            // Check for duplicate samples within this population
            let mut seen = HashSet::new();
            for sample in &samples {
                if !seen.insert(sample) {
                    bail!("Duplicate sample '{}' in population '{}'", sample, pop_name);
                }
            }

            // Build lookup for each sample
            for (sample_idx, sample) in samples.iter().enumerate() {
                if sample_lookup
                    .insert(sample.clone(), (pop_idx, sample_idx))
                    .is_some()
                {
                    bail!("Sample '{}' appears in multiple populations", sample);
                }
            }

            populations.insert(pop_name.clone(), Population::new(pop_name, samples));
        }

        Ok(Self {
            populations,
            sample_lookup,
        })
    }

    /// Validate that this population map matches the given samples exactly
    pub fn validate_exact_match(&self, available_samples: &[String]) -> Result<()> {
        let available: HashSet<_> = available_samples.iter().collect();
        let required: HashSet<_> = self.sample_lookup.keys().collect();

        // Check for samples in pop file but not in data
        let missing: Vec<_> = required
            .difference(&available)
            .map(|s| s.as_str())
            .collect();
        if !missing.is_empty() {
            bail!(
                "Samples in population file not found in data: {}",
                missing.join(", ")
            );
        }

        // Check for samples in data but not in pop file
        let extra: Vec<_> = available
            .difference(&required)
            .map(|s| s.as_str())
            .collect();
        if !extra.is_empty() {
            bail!(
                "Samples in data not found in population file: {}",
                extra.join(", ")
            );
        }

        Ok(())
    }

    /// Get (population_index, sample_index_within_population) for a sample
    pub fn lookup(&self, sample: &str) -> Option<(usize, usize)> {
        self.sample_lookup.get(sample).copied()
    }

    /// Get population by name
    pub fn get_population(&self, name: &str) -> Option<&Population> {
        self.populations.get(name)
    }

    /// Iterate over all populations
    pub fn populations(&self) -> impl Iterator<Item = &Population> {
        self.populations.values()
    }

    /// Get population by index
    pub fn get_population_by_index(&self, idx: usize) -> Option<&Population> {
        self.populations.get_index(idx).map(|(_, pop)| pop)
    }

    /// Number of populations
    pub fn num_populations(&self) -> usize {
        self.populations.len()
    }

    /// All sample names (for validation)
    pub fn all_samples(&self) -> Vec<&str> {
        self.sample_lookup.keys().map(|s| s.as_str()).collect()
    }

    /// Sample counts per population
    pub fn sample_counts_per_population(&self) -> Vec<usize> {
        self.populations.values().map(|p| p.num_samples()).collect()
    }
}
