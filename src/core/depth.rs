use crate::core::contig::{validate_contig_consistency, ContigChunk, ContigSet};
use crate::core::depth::d4::D4Reader;
use color_eyre::eyre::eyre;
use color_eyre::Help;
use color_eyre::Result;
use log::info;
use rayon::prelude::*;
use std::path::{Path, PathBuf};
use indicatif::{ProgressBar, ProgressStyle};
use crate::core::zarr::DepthArrays;
pub mod array;
pub mod d4;
pub mod gvcf;

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

    
}

pub fn open_depth_source(path: impl AsRef<Path>, min_gq: Option<isize>) -> Result<Box<dyn DepthSource>> {
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
            Ok(Box::new(gvcf::GvcfReader::new(path, &sample_name, min_gq)?))
        }
    }
}
pub fn open_depth_sources(paths: Vec<PathBuf>, min_gq: Option<isize>) -> Result<Vec<Box<dyn DepthSource>>> {
    let mut all_sources = Vec::new();

    for path in paths {
        // Auto-detect and expand multisample files
        if is_multisample_d4(&path)? {
            all_sources.extend(D4Reader::from_multisample_file(&path)?);
        } else {
            all_sources.push(open_depth_source(&path, min_gq)?);
        }
    }

    Ok(all_sources)
}

fn is_multisample_d4(path: &Path) -> Result<bool> {
    let path_str = path.to_string_lossy();

    // Only uncompressed .d4 files can be multisample
    if !matches!(DepthFileType::from_path(&path_str), Some(DepthFileType::D4)) {
        return Ok(false);
    }

    let tracks = D4Reader::list_tracks(path)?;
    Ok(tracks.len() > 1)
}

pub enum DepthInput {
    Files(Vec<PathBuf>),
    Zarr(DepthArrays),
}

enum DepthSources {
    IndividualFiles(Vec<PathBuf>),
    MultisampleD4 {
        file: PathBuf,
        tracks: Vec<(String, String)>, // (track_name, sample_name)
    },
}

/// Processor for reading depth data in parallel chunks
pub struct DepthProcessor {
    sources: DepthSources,
    sample_names: Vec<String>,
    reference_contigs: ContigSet,
    min_gq: Option<isize>
}

impl DepthProcessor {
    pub fn from_paths(paths: Vec<PathBuf>, min_gq: Option<isize>) -> Result<Self> {
        info!("Processing {} depth files: {:?}", paths.len(), &paths);
        let (sources, sample_names, reference_contigs) =
            if paths.len() == 1 && is_multisample_d4(&paths[0])? {
                setup_multisample(&paths[0])?
            } else {
                setup_individual_files(paths, min_gq)?
            };
        Ok(Self {
            sources,
            sample_names,
            reference_contigs,
            min_gq
        })
    }

    pub fn sample_names(&self) -> &[String] {
        &self.sample_names
    }

    pub fn reference_contigs(&self) -> &ContigSet {
        &self.reference_contigs
    }

    pub fn process_chunks<F, T>(&self, chunk_size: u64, process_chunk: F) -> Result<Vec<T>>
where
    F: Fn(&ContigChunk, Vec<Vec<u32>>, &[String]) -> Result<T> + Sync,
    T: Send,
{
    let chunks = self.reference_contigs.to_chunks(chunk_size);
    
    let pb = ProgressBar::new(chunks.len() as u64);
    pb.set_style(
        ProgressStyle::default_bar()
            .template("{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {pos}/{len} chunks ({eta})")
            .unwrap()
            .progress_chars("#>-")
    );

    let results: Result<Vec<T>> = chunks
        .par_iter()
        .map(|chunk| {
            let all_sample_depths = self.read_chunk_depths(chunk)?;
            let result = process_chunk(chunk, all_sample_depths, &self.sample_names)?;
            pb.inc(1);
            Ok(result)
        })
        .collect();

    pb.finish_with_message("done");
    results
}

    fn read_chunk_depths(&self, chunk: &ContigChunk) -> Result<Vec<Vec<u32>>> {
        match &self.sources {
            DepthSources::IndividualFiles(paths) => paths
                .iter()
                .map(|path| {
                    let mut reader = open_depth_source(path, self.min_gq)?;
                    reader.read_depths(&chunk.contig_name, chunk.start, chunk.end)
                })
                .collect(),
            DepthSources::MultisampleD4 { file, tracks } => tracks
                .iter()
                .map(|(track_name, sample_name)| {
                    let mut reader = D4Reader::new(file, sample_name, Some(track_name))?;
                    reader.read_depths(&chunk.contig_name, chunk.start, chunk.end)
                })
                .collect(),
        }
    }
}

fn setup_multisample(path: &Path) -> Result<(DepthSources, Vec<String>, ContigSet)> {
    let track_paths = D4Reader::list_tracks(path)?;

    let tracks: Vec<(String, String)> = track_paths
        .into_iter()
        .map(|track_path| {
            let track_name = track_path
                .to_str()
                .ok_or_else(|| eyre!("Invalid track path: {:?}", track_path))?
                .to_string();
            let sample_name = track_name.trim_start_matches('/').to_string();
            Ok((track_name, sample_name))
        })
        .collect::<Result<_>>()?;

    let sample_names: Vec<String> = tracks.iter().map(|(_, name)| name.clone()).collect();

    let first_reader = D4Reader::new(path, &tracks[0].1, Some(&tracks[0].0))?;
    let reference_contigs = first_reader.contigs();

    Ok((
        DepthSources::MultisampleD4 {
            file: path.to_path_buf(),
            tracks,
        },
        sample_names,
        reference_contigs,
    ))
}

fn setup_individual_files(paths: Vec<PathBuf>, min_gq: Option<isize>) -> Result<(DepthSources, Vec<String>, ContigSet)> {
    let readers: Vec<_> = paths
        .iter()
        .map(|path| open_depth_source(path, min_gq))
        .collect::<Result<_>>()?;

    let sample_names: Vec<String> = readers
        .iter()
        .map(|r| r.sample_name().to_string())
        .collect();

    let contig_sets: Vec<_> = paths
        .iter()
        .zip(readers.iter())
        .map(|(path, reader)| (path, reader.contigs()))
        .collect();

    let reference_contigs = validate_contig_consistency(contig_sets)?;

    Ok((
        DepthSources::IndividualFiles(paths),
        sample_names,
        reference_contigs,
    ))
}

#[cfg(test)]
mod tests {
    use super::*;
    use rstest::rstest;

    const TEST_GVCF: &str = "tests/data/depth/all20/sample1.g.vcf.gz";
    const TEST_D4: &str = "tests/data/depth/all20/sample1.d4";
    const TEST_D4_GZ: &str = "tests/data/depth/all20/sample1.d4.gz";
    const TEST_MULTISAMPLE_D4: &str = "tests/data/depth/all20/merged_s1s2.d4";
    const SAMPLE_NAME: &str = "sample1";
    const EXPECTED_DEPTH: u32 = 20;

    #[rstest]
    #[case(TEST_GVCF, "sample1")]
    #[case(TEST_D4, "sample1")]
    #[case(TEST_D4_GZ, "sample1")]
    fn test_open_depth_source(#[case] path: &str, #[case] expected_sample: &str) {
        let result = open_depth_source(path, None);
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
        let result = open_depth_source("test.bam", None);
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
        let reader = open_depth_source(path, None).unwrap();
        assert_eq!(reader.sample_name(), SAMPLE_NAME);
    }

    #[rstest]
    #[case(TEST_GVCF)]
    #[case(TEST_D4)]
    #[case(TEST_D4_GZ)]
    fn test_contigs(#[case] path: &str) {
        let reader = open_depth_source(path, None).unwrap();
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
        let mut reader = open_depth_source(path, None).unwrap();

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
        let mut reader = open_depth_source(path, None).unwrap();

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
        let mut reader = open_depth_source(path, None).unwrap();

        for chrom in ["chr1", "chr2", "chr3", "chr4"] {
            let depths = reader.read_depths(chrom, 0, 1000).unwrap();
            assert_eq!(depths.len(), 1000);

            // All positions should have depth 20
            let sum: u32 = depths.iter().sum();
            assert_eq!(sum, 1000 * EXPECTED_DEPTH);
        }
    }

    #[test]
    fn test_processor_individual_files() {
        let paths = vec![PathBuf::from(TEST_D4), PathBuf::from(TEST_D4_GZ)];
        let processor = DepthProcessor::from_paths(paths, None).unwrap();

        assert_eq!(processor.sample_names(), &["sample1", "sample1"]);
        assert_eq!(processor.reference_contigs().len(), 4);
    }

    #[test]
    fn test_processor_multisample() {
        let paths = vec![PathBuf::from(TEST_MULTISAMPLE_D4)];
        let processor = DepthProcessor::from_paths(paths, None).unwrap();

        assert_eq!(processor.sample_names(), &["sample1", "sample2"]);
        assert_eq!(processor.reference_contigs().len(), 4);
    }

    #[test]
    fn test_processor_process_chunks() {
        let paths = vec![PathBuf::from(TEST_D4)];
        let processor = DepthProcessor::from_paths(paths, None).unwrap();

        let results: Vec<()> = processor
            .process_chunks(100, |chunk, depths, sample_names| {
                assert_eq!(depths.len(), 1);
                assert_eq!(sample_names, &["sample1"]);
                assert_eq!(depths[0].len(), (chunk.end - chunk.start) as usize);
                assert!(depths[0].iter().all(|&d| d == EXPECTED_DEPTH));
                Ok(())
            })
            .unwrap();

        // 4 chromosomes * 1000bp / 100bp chunks = 40 chunks
        assert_eq!(results.len(), 40);
    }

    #[test]
    fn test_processor_multisample_chunks() {
        let paths = vec![PathBuf::from(TEST_MULTISAMPLE_D4)];
        let processor = DepthProcessor::from_paths(paths, None).unwrap();

        let results: Vec<()> = processor
            .process_chunks(100, |chunk, depths, sample_names| {
                assert_eq!(depths.len(), 2);
                assert_eq!(sample_names, &["sample1", "sample2"]);
                assert_eq!(depths[0].len(), (chunk.end - chunk.start) as usize);
                assert!(depths[0].iter().all(|&d| d == EXPECTED_DEPTH));
                Ok(())
            })
            .unwrap();
        assert_eq!(results.len(), 40);
    }
}
