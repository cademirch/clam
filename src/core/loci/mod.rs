use crate::core::depth::{array::MultisampleDepthArray, open_depth_source, DepthSource};
use crate::core::{
    contig::validate_contig_consistency, population::PopulationMap, zarr::CallableArrays,
};
use color_eyre::Result;
use rayon::prelude::*;
use std::path::PathBuf;

pub struct ThresholdConfig {
    pub min_depth: f64,
    pub max_depth: f64,
    pub min_proportion: f64,
    pub mean_depth_range: (f64, f64),
}

pub fn process_individual_files(
    depth_files: Vec<PathBuf>,
    output_path: PathBuf,
    pop_map: PopulationMap,
    thresholds: ThresholdConfig,
    chunk_size: u64,
) -> Result<CallableArrays> {
    let readers: Vec<Box<dyn DepthSource>> = depth_files
        .iter()
        .map(|path| open_depth_source(path))
        .collect::<Result<_>>()?;

    let contig_sets: Vec<_> = depth_files
        .iter()
        .zip(readers.iter())
        .map(|(path, reader)| (path, reader.contigs()))
        .collect();
    let reference_contigs = validate_contig_consistency(contig_sets)?;

    let sample_names: Vec<String> = readers
        .iter()
        .map(|r| r.sample_name().to_string())
        .collect();
    pop_map.validate_exact_match(&sample_names)?;

    let population_names: Vec<String> = pop_map.populations().map(|s| s.to_string()).collect();
    let output_zarr = CallableArrays::create_new(
        output_path,
        reference_contigs.clone(),
        population_names,
        chunk_size,
    )?;

    let chunks = reference_contigs.to_chunks(chunk_size);

    chunks.par_iter().try_for_each(|chunk| {
        let all_sample_depths: Vec<Vec<u32>> = depth_files
            .par_iter()
            .map(|path| {
                let mut reader = open_depth_source(path)?;
                reader.read_depths(&chunk.contig_name, chunk.start, chunk.end)
            })
            .collect::<Result<_>>()?;

        let depth_array =
            MultisampleDepthArray::from_samples(all_sample_depths, sample_names.clone(), &pop_map)?;

        let callable_mask =
            depth_array.apply_sample_thresholds(thresholds.min_depth, thresholds.max_depth);

        let callable_counts = depth_array.count_callable_per_population(
            &callable_mask,
            thresholds.min_proportion,
            thresholds.mean_depth_range,
        )?;

        output_zarr.write_chunk(&chunk.contig_name, chunk.chunk_idx, callable_counts)?;

        Ok::<_, color_eyre::Report>(())
    })?;

    Ok(output_zarr)
}

#[cfg(test)]
mod tests {
    use super::*;
    use rstest::rstest;
    use tempfile::TempDir;

    #[rstest]
    fn test_one_file_callable(
        #[values(
            "tests/data/depth/all20/sample1.g.vcf.gz",
            "tests/data/depth/all20/sample1.d4",
            "tests/data/depth/all20/sample1.d4.gz"
        )]
        path: &str,
        #[values((0.0, 1), (21.0, 0))] threshold_expected: (f64, u8),
    ) {
        let (min_depth, expected) = threshold_expected;
        let dir = TempDir::new().unwrap();
        let output_path = dir.path().join("callable.zarr");

        let depth_files = vec![PathBuf::from(path)];
        let pop_map = PopulationMap::default_from_samples(vec!["sample1".to_string()]);

        let thresholds = ThresholdConfig {
            min_depth,
            max_depth: f64::INFINITY,
            min_proportion: 0.0,
            mean_depth_range: (0.0, f64::INFINITY),
        };

        let result =
            process_individual_files(depth_files, output_path.clone(), pop_map, thresholds, 100)
                .unwrap();

        assert!(
            output_path.exists(),
            "Output zarr directory should be created"
        );

        for (chrom_name, chrom_length) in result.contigs().iter() {
            let num_chunks = (chrom_length as u64 + 99) / 100;
            for chunk_idx in 0..num_chunks {
                let chunk_data = result.read_chunk(chrom_name, chunk_idx).unwrap();
                assert!(chunk_data.iter().all(|&x| x == expected));
            }
        }
    }
}
