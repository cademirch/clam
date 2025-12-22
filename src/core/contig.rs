use std::path::PathBuf;
use color_eyre::{
    eyre::{bail, WrapErr},
    Result,
};
use indexmap::IndexMap;
use std::collections::HashSet;
/// Information about a reference contig
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Contig {
    pub name: String,
    pub length: usize,
}

impl Contig {
    pub fn new(name: String, length: usize) -> Self {
        Self { name, length }
    }
    
}

#[derive(Debug, Clone)]
pub struct ContigChunk {
    pub contig_name: String,
    pub chunk_idx: u64,
    pub start: u32,
    pub end: u32,
}



#[derive(Debug, Clone)]
pub struct ContigSet {
    contigs: IndexMap<String, usize>, // name -> length, ordered
}

impl ContigSet {
    /// Create from a list of contigs
    pub fn new(contigs: Vec<Contig>) -> Self {
        let contigs = contigs.into_iter().map(|c| (c.name, c.length)).collect();
        Self { contigs }
    }
    
    /// Get contig length by name
    pub fn get_length(&self, name: &str) -> Option<usize> {
        self.contigs.get(name).copied()
    }

    /// Check if a contig exists
    pub fn contains(&self, name: &str) -> bool {
        self.contigs.contains_key(name)
    }

    /// Iterate over contigs in order
    pub fn iter(&self) -> impl Iterator<Item = (&str, usize)> {
        self.contigs
            .iter()
            .map(|(name, &length)| (name.as_str(), length))
    }

    /// Number of contigs
    pub fn len(&self) -> usize {
        self.contigs.len()
    }

    /// Contig names
    pub fn names(&self) -> Vec<&str> {
        self.contigs.keys().map(|s| s.as_str()).collect()
    }

    /// Validate exact match with another contig set
    pub fn validate_matches(&self, other: &ContigSet) -> Result<()> {
        // Check same number of contigs
        if self.contigs.len() != other.contigs.len() {
            bail!(
                "Different number of contigs: {} vs {}",
                self.contigs.len(),
                other.contigs.len()
            );
        }

        // Check each contig matches
        for (name, &length) in &self.contigs {
            match other.get_length(name) {
                None => bail!("Contig '{}' not found in comparison set", name),
                Some(other_length) if other_length != length => {
                    bail!(
                        "Contig '{}' has mismatched lengths: {} vs {}",
                        name,
                        length,
                        other_length
                    );
                }
                Some(_) => {} // OK
            }
        }

        Ok(())
    }

    pub fn to_chunks(&self, chunk_size: u64) -> Vec<ContigChunk> {
        let chunks: Vec<ContigChunk> = self.iter()
            .flat_map(|(contig, length)| {
                let num_chunks = (length + chunk_size as usize - 1) / chunk_size as usize;
                (0..num_chunks)
                    .map(|chunk_idx| {
                        let start = (chunk_idx * chunk_size as usize) as u32;
                        let end = ((chunk_idx + 1) * chunk_size as usize).min(length) as u32;
                        ContigChunk {
                            contig_name: contig.to_string(),
                            chunk_idx: chunk_idx.try_into().unwrap(),
                            start,
                            end,
                        }
                    })
                    .collect::<Vec<_>>()
            })
            .collect();
        
        chunks
    }
    /// Filter contigs based on include/exclude lists
    pub fn filter(
        self,
        include: Option<&HashSet<String>>,
        exclude: Option<&HashSet<String>>,
    ) -> Self {
        let filtered = self.contigs
            .into_iter()
            .filter(|(name, _)| {
                
                if let Some(include_set) = include {
                    if !include_set.contains(name.as_str()) {
                        return false;
                    }
                }
                
                
                if let Some(exclude_set) = exclude {
                    if exclude_set.contains(name.as_str()) {
                        return false;
                    }
                }
                
                true
            })
            .collect();
        
        Self { contigs: filtered }
    }
    
    /// Check if empty after filtering
    pub fn is_empty(&self) -> bool {
        self.contigs.is_empty()
    }
}

pub fn validate_contig_consistency(contig_sets: Vec<(&PathBuf, ContigSet)>) -> Result<ContigSet> {
    if contig_sets.is_empty() {
        bail!("No contig sets provided for validation");
    }

    // Use the first contig set as reference
    let (first_path, reference) = &contig_sets[0];

    // Validate all other contig sets against the reference
    for (path, contig_set) in contig_sets.iter().skip(1) {
        reference.validate_matches(contig_set).wrap_err_with(|| {
            format!(
                "Contig mismatch between '{}' and '{}'",
                first_path.display(),
                path.display()
            )
        })?;
    }

    Ok(reference.clone())
}
