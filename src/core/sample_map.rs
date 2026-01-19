use color_eyre::eyre::eyre;
use color_eyre::Result;
use std::collections::HashMap;
use std::path::{Path, PathBuf};

/// Maps filenames to custom sample names
#[derive(Debug, Clone)]
pub struct SampleMap {
    /// Maps filename (as string) -> sample_name
    filename_to_sample: HashMap<String, String>,
}

impl SampleMap {
    /// Create from a tab-separated file: filename\tsample_name (no header)
    pub fn from_file(path: impl AsRef<Path>) -> Result<Self> {
        let path = path.as_ref();
        let content = std::fs::read_to_string(path)?;

        let mut filename_to_sample = HashMap::new();
        let mut sample_names_seen = HashMap::new();

        for (line_num, line) in content.lines().enumerate() {
            let line = line.trim();
            if line.is_empty() {
                continue;
            }

            let parts: Vec<&str> = line.split('\t').collect();
            if parts.len() != 2 {
                return Err(eyre!(
                    "Invalid line {} in sample file '{}': expected 2 tab-separated columns, got {}",
                    line_num + 1,
                    path.display(),
                    parts.len()
                ));
            }

            let filename = parts[0].to_string();
            let sample_name = parts[1].to_string();

            // Check for duplicate filenames
            if filename_to_sample.contains_key(&filename) {
                return Err(eyre!(
                    "Duplicate filename '{}' in sample file '{}'",
                    filename,
                    path.display()
                ));
            }

            // Check for duplicate sample names
            if let Some(existing_file) = sample_names_seen.get(&sample_name) {
                return Err(eyre!(
                    "Duplicate sample name '{}' in sample file '{}' (used by both '{}' and '{}')",
                    sample_name,
                    path.display(),
                    existing_file,
                    filename
                ));
            }

            sample_names_seen.insert(sample_name.clone(), filename.clone());
            filename_to_sample.insert(filename, sample_name);
        }

        if filename_to_sample.is_empty() {
            return Err(eyre!(
                "Sample file '{}' is empty or contains no valid entries",
                path.display()
            ));
        }

        Ok(Self { filename_to_sample })
    }

    /// Get sample name for a filename, returning None if not in map
    pub fn get_sample_name(&self, filename: &str) -> Option<&str> {
        self.filename_to_sample.get(filename).map(|s| s.as_str())
    }

    /// Validate that all input files are covered by the sample map
    pub fn validate_coverage(&self, input_files: &[PathBuf]) -> Result<()> {
        let missing: Vec<_> = input_files
            .iter()
            .filter_map(|path| {
                let filename = path.file_name()?.to_str()?;
                if self.filename_to_sample.contains_key(filename) {
                    None
                } else {
                    Some(filename.to_string())
                }
            })
            .collect();

        if !missing.is_empty() {
            return Err(eyre!(
                "The following input files are not listed in the sample file: {}",
                missing.join(", ")
            ));
        }

        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;
    use tempfile::NamedTempFile;

    #[test]
    fn test_parse_sample_file() {
        let mut file = NamedTempFile::new().unwrap();
        writeln!(file, "L.sample1.d4.gz\tL.sample1").unwrap();
        writeln!(file, "L.sample2.d4.gz\tL.sample2").unwrap();

        let map = SampleMap::from_file(file.path()).unwrap();

        assert_eq!(map.get_sample_name("L.sample1.d4.gz"), Some("L.sample1"));
        assert_eq!(map.get_sample_name("L.sample2.d4.gz"), Some("L.sample2"));
        assert_eq!(map.get_sample_name("unknown.d4"), None);
    }

    #[test]
    fn test_validate_coverage_success() {
        let mut file = NamedTempFile::new().unwrap();
        writeln!(file, "sample1.d4\tsample1").unwrap();
        writeln!(file, "sample2.d4\tsample2").unwrap();

        let map = SampleMap::from_file(file.path()).unwrap();
        let files = vec![
            PathBuf::from("/path/to/sample1.d4"),
            PathBuf::from("/other/sample2.d4"),
        ];

        assert!(map.validate_coverage(&files).is_ok());
    }

    #[test]
    fn test_validate_coverage_missing_file() {
        let mut file = NamedTempFile::new().unwrap();
        writeln!(file, "sample1.d4\tsample1").unwrap();

        let map = SampleMap::from_file(file.path()).unwrap();
        let files = vec![
            PathBuf::from("sample1.d4"),
            PathBuf::from("sample2.d4"), // Not in map
        ];

        let result = map.validate_coverage(&files);
        assert!(result.is_err());
        assert!(result.unwrap_err().to_string().contains("sample2.d4"));
    }

    #[test]
    fn test_duplicate_filename_error() {
        let mut file = NamedTempFile::new().unwrap();
        writeln!(file, "sample1.d4\tsample1").unwrap();
        writeln!(file, "sample1.d4\tsample1_dup").unwrap(); // Duplicate filename

        let result = SampleMap::from_file(file.path());
        assert!(result.is_err());
        assert!(result.unwrap_err().to_string().contains("Duplicate filename"));
    }

    #[test]
    fn test_duplicate_sample_name_error() {
        let mut file = NamedTempFile::new().unwrap();
        writeln!(file, "sample1.d4\tshared_name").unwrap();
        writeln!(file, "sample2.d4\tshared_name").unwrap(); // Duplicate sample name

        let result = SampleMap::from_file(file.path());
        assert!(result.is_err());
        assert!(result
            .unwrap_err()
            .to_string()
            .contains("Duplicate sample name"));
    }

    #[test]
    fn test_empty_file_error() {
        let file = NamedTempFile::new().unwrap();
        // File is empty

        let result = SampleMap::from_file(file.path());
        assert!(result.is_err());
        assert!(result.unwrap_err().to_string().contains("empty"));
    }

    #[test]
    fn test_invalid_line_format() {
        let mut file = NamedTempFile::new().unwrap();
        writeln!(file, "only_one_column").unwrap();

        let result = SampleMap::from_file(file.path());
        assert!(result.is_err());
        assert!(result.unwrap_err().to_string().contains("2 tab-separated"));
    }

    #[test]
    fn test_skips_empty_lines() {
        let mut file = NamedTempFile::new().unwrap();
        writeln!(file, "sample1.d4\tsample1").unwrap();
        writeln!(file).unwrap(); // Empty line
        writeln!(file, "   ").unwrap(); // Whitespace only
        writeln!(file, "sample2.d4\tsample2").unwrap();

        let map = SampleMap::from_file(file.path()).unwrap();
        assert_eq!(map.get_sample_name("sample1.d4"), Some("sample1"));
        assert_eq!(map.get_sample_name("sample2.d4"), Some("sample2"));
    }
}
