use crate::core::contig::ContigSet;
use color_eyre::eyre::eyre;
use color_eyre::Help;
use color_eyre::Result;
use std::path::Path;

/// Trait for reading depth information from files
pub trait DepthSource {
    /// Read raw depth values for a region
    fn read_depths(&mut self, chrom: &str, start: u32, end: u32) -> Result<Vec<u32>>;

    /// Get contigs available in this file
    fn contigs(&self) -> ContigSet;

    /// Get the sample name for this reader
    fn sample_name(&self) -> &str;
}

enum DepthFileType {
    D4,
    D4Gz,
    GvcfGz,
}

impl DepthFileType {
    fn from_path(path: &str) -> Option<Self> {
        if path.ends_with(".d4.gz") {
            Some(Self::D4Gz)
        } else if path.ends_with(".d4") {
            Some(Self::D4)
        } else if path.ends_with(".g.vcf.gz") || path.ends_with(".gvcf.gz") {
            Some(Self::GvcfGz)
        } else {
            None
        }
    }

    fn suffix(&self) -> &str {
        match self {
            Self::D4 => ".d4",
            Self::D4Gz => ".d4.gz",
            Self::GvcfGz => ".g.vcf.gz",
        }
    }
}

pub fn open_depth_source(path: impl AsRef<Path>) -> Result<Box<dyn DepthSource>> {
    let path = path.as_ref();
    let path_str = path.to_string_lossy();

    let file_type = DepthFileType::from_path(&path_str).ok_or_else(|| {
        eyre!(
            "Unsupported depth file format: {}. Supported formats: .d4, .d4.gz, .g.vcf.gz, .gvcf.gz",
            path_str
        )
    })?;

    let sample_name = path
        .file_prefix()
        .and_then(|s| s.to_str())
        .ok_or_else(|| eyre!("Could not extract sample name from path: {}", path_str))
        .suggestion(
            "Sample names are derived from file prefixes (e.g., 'sample' from 'sample.d4.gz')",
        )?
        .to_string();
    match file_type {
        DepthFileType::D4Gz => {
            // BGZF-compressed D4 file
            Ok(Box::new(d4::BgzfD4Reader::new(path, &sample_name, None)?))
        }
        DepthFileType::D4 => {
            // Uncompressed D4 file
            Ok(Box::new(d4::D4Reader::new(path, &sample_name, None)?))
        }
        DepthFileType::GvcfGz => {
            // GVCF file
            Ok(Box::new(gvcf::GvcfReader::new(path, &sample_name)?))
        }
    }
}

pub mod array;
mod d4;
mod gvcf;

#[cfg(test)]
mod tests {
    use super::*;
    use rstest::rstest;

    const TEST_GVCF: &str = "tests/data/depth/all20/sample1.g.vcf.gz";
    const TEST_D4: &str = "tests/data/depth/all20/sample1.d4";
    const TEST_D4_GZ: &str = "tests/data/depth/all20/sample1.d4.gz";
    const SAMPLE_NAME: &str = "sample1";
    const EXPECTED_DEPTH: u32 = 20;

    #[rstest]
    #[case(TEST_GVCF, "sample1")]
    #[case(TEST_D4, "sample1")]
    #[case(TEST_D4_GZ, "sample1")]
    fn test_open_depth_source(#[case] path: &str, #[case] expected_sample: &str) {
        let result = open_depth_source(path);
        assert!(
            result.is_ok(),
            "Failed to open {}: {:?}",
            path,
            result.err()
        );

        let reader = result.unwrap();
        assert_eq!(reader.sample_name(), expected_sample);
    }

    #[test]
    fn test_unsupported_format() {
        let result = open_depth_source("test.bam");
        assert!(result.is_err());

        let err = result.err().unwrap();
        let err_msg = err.to_string();
        assert!(err_msg.contains("Unsupported depth file format"));
    }

    #[test]
    fn test_file_type_detection() {
        assert!(matches!(
            DepthFileType::from_path("sample.d4"),
            Some(DepthFileType::D4)
        ));
        assert!(matches!(
            DepthFileType::from_path("sample.d4.gz"),
            Some(DepthFileType::D4Gz)
        ));
        assert!(matches!(
            DepthFileType::from_path("sample.g.vcf.gz"),
            Some(DepthFileType::GvcfGz)
        ));
        assert!(matches!(
            DepthFileType::from_path("sample.gvcf.gz"),
            Some(DepthFileType::GvcfGz)
        ));
        assert!(DepthFileType::from_path("sample.bam").is_none());
    }

    // Generic DepthSource tests using rstest matrix
    #[rstest]
    #[case(TEST_GVCF)]
    #[case(TEST_D4)]
    #[case(TEST_D4_GZ)]
    fn test_sample_name(#[case] path: &str) {
        let reader = open_depth_source(path).unwrap();
        assert_eq!(reader.sample_name(), SAMPLE_NAME);
    }

    #[rstest]
    #[case(TEST_GVCF)]
    #[case(TEST_D4)]
    #[case(TEST_D4_GZ)]
    fn test_contigs(#[case] path: &str) {
        let reader = open_depth_source(path).unwrap();
        let contigs = reader.contigs();

        // Should have 4 chromosomes
        assert_eq!(contigs.len(), 4);

        // Each should be 1000 bp
        for contig in contigs.iter() {
            assert_eq!(contig.1, 1000);
        }
    }

    #[rstest]
    #[case(TEST_GVCF)]
    #[case(TEST_D4)]
    #[case(TEST_D4_GZ)]
    fn test_read_depths_full_region(#[case] path: &str) {
        let mut reader = open_depth_source(path).unwrap();

        // Read first 100 positions of chr1
        let depths = reader.read_depths("chr1", 0, 100).unwrap();

        assert_eq!(depths.len(), 100);
        for depth in &depths {
            assert_eq!(*depth, EXPECTED_DEPTH);
        }
    }

    #[rstest]
    #[case(TEST_GVCF)]
    #[case(TEST_D4)]
    #[case(TEST_D4_GZ)]
    fn test_read_depths_partial_region(#[case] path: &str) {
        let mut reader = open_depth_source(path).unwrap();

        // Read a smaller region
        let depths = reader.read_depths("chr2", 100, 200).unwrap();

        assert_eq!(depths.len(), 100);
        for depth in &depths {
            assert_eq!(*depth, EXPECTED_DEPTH);
        }
    }

    #[rstest]
    #[case(TEST_GVCF)]
    #[case(TEST_D4)]
    #[case(TEST_D4_GZ)]
    fn test_read_depths_all_chromosomes(#[case] path: &str) {
        let mut reader = open_depth_source(path).unwrap();

        for chrom in ["chr1", "chr2", "chr3", "chr4"] {
            let depths = reader.read_depths(chrom, 0, 1000).unwrap();
            assert_eq!(depths.len(), 1000);

            // All positions should have depth 20
            let sum: u32 = depths.iter().sum();
            assert_eq!(sum, 1000 * EXPECTED_DEPTH);
        }
    }
}
