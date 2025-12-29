use crate::core::depth::array::stack_depths;
use crate::core::depth::DepthProcessor;
use crate::core::zarr::{DepthArrays, is_zarr_path};
use color_eyre::Result;
use color_eyre::Help;
use color_eyre::eyre::eyre;
use std::path::PathBuf;


pub fn run_collect(depth_files: Vec<PathBuf>, output_path: PathBuf, chunk_size: u64, min_gq: Option<isize>) -> Result<()> {
    if is_zarr_path(&output_path) {
        return Err(eyre!(
            "Output zarr path: {} already exists",
            output_path.display()
        ))
        .suggestion("Remove the existing directory or choose a different output path");
    }
    
    let processor = DepthProcessor::from_paths(depth_files, min_gq)?;

    let output_zarr = DepthArrays::create_new(
        output_path,
        processor.reference_contigs().clone(),
        processor.sample_names().to_vec(),
        chunk_size,
    )?;

    processor.process_chunks(chunk_size, |chunk, depths, sample_names| {
        let depth_array = stack_depths(depths, sample_names.to_vec())?;
        output_zarr.write_chunk(&chunk.contig_name, chunk.chunk_idx, depth_array, None)?;
        Ok(())
    })?;

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use rstest::rstest;
    use tempfile::TempDir;

    #[rstest]
    fn test_collect_one_file(
        #[values(
            "tests/data/depth/all20/sample1.g.vcf.gz",
            "tests/data/depth/all20/sample1.d4",
            "tests/data/depth/all20/sample1.d4.gz"
        )]
        path: &str,
    ) {
        let dir = TempDir::new().unwrap();
        let output_path = dir.path().join("depths.zarr");

        let depth_files = vec![PathBuf::from(path)];

        run_collect(depth_files, output_path.clone(), 100, None).unwrap();

        assert!(output_path.exists());

        let result = DepthArrays::open(&output_path).unwrap();
        assert_eq!(result.column_names(), &["sample1"]);

        for (chrom_name, chrom_length) in result.contigs().iter() {
            let num_chunks = (chrom_length as u64 + 99) / 100;
            for chunk_idx in 0..num_chunks {
                let chunk_data = result.read_chunk(chrom_name, chunk_idx, None).unwrap();

                let expected_rows = std::cmp::min(100, chrom_length - (chunk_idx * 100) as usize);
                assert_eq!(chunk_data.shape(), &[expected_rows, 1]);

                assert!(chunk_data.iter().all(|&x| x == 20));
            }
        }
    }

    #[rstest]
    fn test_collect_two_samples(
        #[values(
            vec!["tests/data/depth/all20/sample1.d4", "tests/data/depth/all20/sample2.d4"],
            vec!["tests/data/depth/all20/merged_s1s2.d4"]
        )]
        paths: Vec<&str>,
    ) {
        let dir = TempDir::new().unwrap();
        let output_path = dir.path().join("depths.zarr");

        let depth_files: Vec<PathBuf> = paths.iter().map(|p| PathBuf::from(p)).collect();

        run_collect(depth_files, output_path.clone(), 100, None).unwrap();

        assert!(output_path.exists());

        let result = DepthArrays::open(&output_path).unwrap();
        assert_eq!(result.column_names(), &["sample1", "sample2"]);

        for (chrom_name, chrom_length) in result.contigs().iter() {
            let num_chunks = (chrom_length as u64 + 99) / 100;
            for chunk_idx in 0..num_chunks {
                let chunk_data = result.read_chunk(chrom_name, chunk_idx, None).unwrap();

                let expected_rows = std::cmp::min(100, chrom_length - (chunk_idx * 100) as usize);
                assert_eq!(chunk_data.shape(), &[expected_rows, 2]);

                assert!(chunk_data.iter().all(|&x| x == 20));
            }
        }
    }
}
