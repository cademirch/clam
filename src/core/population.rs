use color_eyre::{
    eyre::{bail, Context, ContextCompat},
    Result,
};
use indexmap::IndexMap;
use log::warn;
use ndarray::Array2;
use serde::{Deserialize, Serialize};
use std::collections::{HashMap, HashSet};
use std::fmt;
use std::fs::File;
use std::path::{Path, PathBuf};

/// A single population with its samples
#[derive(Debug, Clone, Serialize, Deserialize)]
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

    pub fn from_populations(pop_data: IndexMap<String, Vec<String>>) -> Result<Self> {
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
    ///
    /// If `force_samples` is true, samples in data but not in population file
    /// will be excluded from analysis (different warning message).
    pub fn validate_exact_match(
        &self,
        available_samples: &[String],
        force_samples: bool,
    ) -> Result<()> {
        let available: HashSet<_> = available_samples.iter().collect();
        let required: HashSet<_> = self.sample_lookup.keys().collect();

        // Check for samples in pop file but not in data
        let missing: Vec<_> = required
            .difference(&available)
            .map(|s| s.as_str())
            .collect();
        if !missing.is_empty() {
            warn!(
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
            let display = if extra.len() <= 5 {
                extra.join(", ")
            } else {
                format!(
                    "{}, ... and {} more",
                    extra[..5].join(", "),
                    extra.len() - 5
                )
            };

            if force_samples {
                warn!(
                    "{} sample(s) will be excluded from analysis (not in population file): {}",
                    extra.len(),
                    display
                );
            } else {
                warn!(
                    "{} samples in data not found in population file: {}",
                    extra.len(),
                    display
                );
            }
        }

        Ok(())
    }

    /// Create a population membership matrix for the given samples
    ///
    /// Returns an Array2 where:
    /// - Rows represent samples (in the order provided)
    /// - Columns represent populations (in the order they were added)
    /// - Values are 1 if sample belongs to that population, 0 otherwise
    ///
    /// # Example
    /// For 3 samples across 2 populations:
    /// ```text
    /// [[1, 0],  // sample 0 in pop 0
    ///  [0, 1],  // sample 1 in pop 1
    ///  [1, 0]]  // sample 2 in pop 0
    /// ```
    pub fn membership_matrix(&self, samples: &[String]) -> Result<Array2<usize>> {
        let num_samples = samples.len();
        let num_pops = self.num_populations();

        let mut data = vec![0usize; num_samples * num_pops];

        for (sample_idx, sample) in samples.iter().enumerate() {
            let (pop_idx, _) = self
                .lookup(sample)
                .wrap_err_with(|| format!("Sample '{}' not found in population map", sample))?;

            // Set the appropriate element to 1 (row-major order)
            data[sample_idx * num_pops + pop_idx] = 1;
        }

        Array2::from_shape_vec((num_samples, num_pops), data)
            .wrap_err("Failed to create membership matrix")
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

    /// Get vector of owned poulations
    pub fn populations_owned(&self) -> Vec<Population> {
        self.populations().cloned().collect()
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

    /// Filter a sample list to only include samples in this population map.
    /// Returns the filtered samples in the same order as the input.
    pub fn filter_samples(&self, samples: &[String]) -> Vec<String> {
        samples
            .iter()
            .filter(|s| self.sample_lookup.contains_key(*s))
            .cloned()
            .collect()
    }

    /// Get indices of samples that are in this population map.
    /// Returns Vec of indices (into the input samples array) that should be kept.
    pub fn filter_sample_indices(&self, samples: &[String]) -> Vec<usize> {
        samples
            .iter()
            .enumerate()
            .filter(|(_, s)| self.sample_lookup.contains_key(*s))
            .map(|(idx, _)| idx)
            .collect()
    }
}

/// Row from samples TSV file (for serde deserialization)
#[derive(Debug, Clone, Deserialize)]
struct SampleRow {
    sample_name: String,
    #[serde(default)]
    file_path: Option<String>,
    #[serde(default)]
    population: Option<String>,
}

/// Entry from samples TSV file
#[derive(Debug, Clone)]
pub struct SampleEntry {
    pub sample_name: String,
    pub file_path: Option<PathBuf>,
    pub population: Option<String>,
}

/// Unified samples configuration from --samples TSV
///
/// TSV format (with header):
/// ```text
/// sample_name    file_path    population
/// sample1        /path/to/sample1.d4.gz    pop_A
/// sample2        /path/to/sample2.d4.gz    pop_B
/// ```
///
/// - `sample_name` column is required
/// - `file_path` column is required when `require_file_path` is true (e.g. for loci/collect), optional otherwise (e.g. for stat)
/// - `population` column is optional; if present, ALL rows must have values; if absent, all samples assigned to "default" population
#[derive(Debug, Clone)]
pub struct SamplesConfig {
    entries: Vec<SampleEntry>,
    has_populations: bool,
}

impl SamplesConfig {
    /// Parse from a TSV file with header row.
    ///
    /// When `require_file_path` is true, the `file_path` column must be present
    /// and non-empty for every row. When false, it is optional (useful for
    /// subcommands like `stat` that only need sample/population info).
    pub fn from_file(path: impl AsRef<Path>, require_file_path: bool) -> Result<Self> {
        let path = path.as_ref();
        let file = File::open(path)
            .wrap_err_with(|| format!("Failed to open samples file: {}", path.display()))?;

        let mut reader = csv::ReaderBuilder::new()
            .delimiter(b'\t')
            .has_headers(true)
            .flexible(true)
            .from_reader(file);

        // Verify required columns exist
        let headers = reader.headers().wrap_err("Failed to read header row")?;
        if !headers.iter().any(|h| h == "sample_name") {
            bail!("Missing required 'sample_name' column in header");
        }
        let has_file_path_column = headers.iter().any(|h| h == "file_path");
        if require_file_path && !has_file_path_column {
            bail!("Missing required 'file_path' column in header");
        }

        let mut entries = Vec::new();
        let mut seen_samples: HashSet<String> = HashSet::new();
        let mut seen_filenames: HashMap<String, String> = HashMap::new();

        for (line_num, result) in reader.deserialize().enumerate() {
            let row: SampleRow = result.wrap_err_with(|| {
                format!("Failed to parse line {} in samples file", line_num + 2)
            })?;

            if row.sample_name.is_empty() {
                bail!("Empty sample_name on line {}", line_num + 2);
            }

            // Normalize empty file_path to None
            let file_path_str = row.file_path.filter(|p| !p.is_empty());

            if require_file_path && file_path_str.is_none() {
                bail!("Empty file_path on line {}", line_num + 2);
            }

            // Check for duplicate sample names
            if seen_samples.contains(&row.sample_name) {
                bail!(
                    "Duplicate sample name '{}' on line {}",
                    row.sample_name,
                    line_num + 2
                );
            }
            seen_samples.insert(row.sample_name.clone());

            // Convert file_path and check for duplicate filenames (only when paths are present)
            let file_path = file_path_str.map(|p| PathBuf::from(&p));
            if let Some(ref fp) = file_path {
                if let Some(filename) = fp.file_name().and_then(|f| f.to_str()) {
                    if let Some(existing) = seen_filenames.get(filename) {
                        bail!(
                            "Duplicate filename '{}' for samples '{}' and '{}'",
                            filename,
                            existing,
                            row.sample_name
                        );
                    }
                    seen_filenames.insert(filename.to_string(), row.sample_name.clone());
                }
            }

            // Handle empty string as None for population
            let population = row.population.filter(|p| !p.is_empty());

            entries.push(SampleEntry {
                sample_name: row.sample_name,
                file_path,
                population,
            });
        }

        if entries.is_empty() {
            bail!("Samples file is empty (no data rows)");
        }

        // Validate consistency: if any entry has population, all must have it
        let has_populations = entries.iter().any(|e| e.population.is_some());
        if has_populations && entries.iter().any(|e| e.population.is_none()) {
            bail!("If population column has values, ALL rows must have population values");
        }

        Ok(Self {
            entries,
            has_populations,
        })
    }

    /// Create from file paths, deriving sample names from filename prefixes.
    /// Used for legacy positional args path.
    pub fn from_paths(paths: Vec<PathBuf>) -> Result<Self> {
        use color_eyre::Help;

        let mut entries = Vec::new();
        let mut seen_samples: HashSet<String> = HashSet::new();
        let mut seen_filenames: HashMap<String, String> = HashMap::new();

        for path in paths {
            let sample_name = path
                .file_prefix()
                .and_then(|s| s.to_str())
                .map(|s| s.to_string())
                .ok_or_else(|| {
                    color_eyre::eyre::eyre!("Could not extract sample name from path: {}", path.display())
                })
                .suggestion(
                    "Sample names are derived from file prefixes (e.g., 'sample' from 'sample.d4.gz')",
                )?;

            if !seen_samples.insert(sample_name.clone()) {
                bail!("Duplicate sample name '{}' derived from paths", sample_name);
            }

            // Check for duplicate filenames
            if let Some(filename) = path.file_name().and_then(|f| f.to_str()) {
                if let Some(existing) = seen_filenames.get(filename) {
                    bail!(
                        "Duplicate filename '{}' for samples '{}' and '{}'",
                        filename,
                        existing,
                        sample_name
                    );
                }
                seen_filenames.insert(filename.to_string(), sample_name.clone());
            }

            entries.push(SampleEntry {
                sample_name,
                file_path: Some(path),
                population: None,
            });
        }

        if entries.is_empty() {
            bail!("No input files provided");
        }

        Ok(Self {
            entries,
            has_populations: false,
        })
    }

    /// Get ordered sample names
    pub fn sample_names(&self) -> Vec<String> {
        self.entries.iter().map(|e| e.sample_name.clone()).collect()
    }

    /// Get file paths. Panics if any entry is missing a file_path.
    /// Only call this when the config was created with `require_file_path = true`
    /// or via `from_paths()`.
    pub fn file_paths(&self) -> Vec<PathBuf> {
        self.entries
            .iter()
            .map(|e| {
                e.file_path.clone().expect(
                    "file_path is None; was SamplesConfig created with require_file_path=true?",
                )
            })
            .collect()
    }

    /// Check if population column is present and populated
    pub fn has_populations(&self) -> bool {
        self.has_populations
    }

    /// Convert to PopulationMap
    /// If no population column, all samples are assigned to "default" population
    pub fn to_population_map(&self) -> PopulationMap {
        let mut pop_data: IndexMap<String, Vec<String>> = IndexMap::new();

        if self.has_populations {
            for entry in &self.entries {
                let pop_name = entry.population.as_ref().unwrap();
                pop_data
                    .entry(pop_name.clone())
                    .or_default()
                    .push(entry.sample_name.clone());
            }
        } else {
            // All samples in "default" population
            let all_samples: Vec<String> =
                self.entries.iter().map(|e| e.sample_name.clone()).collect();
            pop_data.insert("default".to_string(), all_samples);
        }

        PopulationMap::from_populations(pop_data).expect("Population map should be valid")
    }

    /// Get number of samples
    pub fn len(&self) -> usize {
        self.entries.len()
    }

    /// Check if empty
    pub fn is_empty(&self) -> bool {
        self.entries.is_empty()
    }
}
