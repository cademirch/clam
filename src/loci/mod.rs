use crate::core::depth::array::{build_pop_membership, MultisampleDepthArray};
use crate::core::depth::DepthProcessor;
use crate::core::population::{PopulationMap, SamplesConfig};
use crate::core::utils::{create_progress_bar, create_spinner};
use crate::core::zarr::{CallableArrays, DepthArrays, SampleMaskArrays};
use color_eyre::Result;
use rayon::iter::{IntoParallelIterator, ParallelIterator};
use std::collections::HashMap;
use std::path::PathBuf;

pub struct ThresholdConfig {
    pub min_depth: f64,
    pub max_depth: f64,
    pub min_proportion: f64,
    pub mean_depth_range: (f64, f64),
    pub per_contig: Option<HashMap<String, (f64, f64)>>, // contig -> (min, max)
}

impl ThresholdConfig {
    pub fn depth_range_for_contig(&self, contig: &str) -> (f64, f64) {
        self.per_contig
            .as_ref()
            .and_then(|m| m.get(contig))
            .copied()
            .unwrap_or((self.min_depth, self.max_depth))
    }
}

/// Run loci calculation from a SamplesConfig
pub fn run_loci(
    config: &SamplesConfig,
    output_path: PathBuf,
    pop_map: Option<PopulationMap>,
    thresholds: ThresholdConfig,
    chunk_size: u64,
    output_per_sample_mask: bool,
    min_gq: Option<isize>,
) -> Result<()> {
    let processor = DepthProcessor::from_samples_config(config, min_gq)?;

    let spinner = create_spinner("Validating population map...");
    //TODO: check that threshold per contig matches depth processor contigs
    let pop_map = pop_map.unwrap_or_else(|| config.to_population_map());
    pop_map.validate_exact_match(processor.sample_names(), false)?;
    spinner.finish_and_clear();

    if output_per_sample_mask {
        process_sample_masks(processor, output_path, thresholds, chunk_size, pop_map)?;
    } else {
        process_population_counts(processor, output_path, pop_map, thresholds, chunk_size)?;
    }

    Ok(())
}

fn process_sample_masks(
    processor: DepthProcessor,
    output_path: PathBuf,
    thresholds: ThresholdConfig,
    chunk_size: u64,
    pop_map: PopulationMap,
) -> Result<SampleMaskArrays> {
    let output_zarr = SampleMaskArrays::create_new(
        output_path,
        processor.reference_contigs().clone(),
        processor.sample_names().to_vec(),
        chunk_size,
        Some(pop_map.populations_owned()),
    )?;

    processor.process_chunks(chunk_size, |chunk, depths, sample_names| {
        let depth_array = MultisampleDepthArray::from_samples(depths, sample_names.to_vec())?;
        let callable_mask =
            depth_array.apply_sample_thresholds(thresholds.min_depth, thresholds.max_depth);
        output_zarr.write_chunk(&chunk.contig_name, chunk.chunk_idx, callable_mask, None)?;
        Ok(())
    })?;

    Ok(output_zarr)
}

fn process_population_counts(
    processor: DepthProcessor,
    output_path: PathBuf,
    pop_map: PopulationMap,
    thresholds: ThresholdConfig,
    chunk_size: u64,
) -> Result<CallableArrays> {
    let population_names: Vec<String> = pop_map.populations().map(|p| p.name.clone()).collect();

    let output_zarr = CallableArrays::create_new(
        output_path,
        processor.reference_contigs().clone(),
        population_names,
        chunk_size,
        Some(pop_map.populations_owned()),
    )?;

    let pop_membership = build_pop_membership(processor.sample_names(), &pop_map);

    processor.process_chunks(chunk_size, |chunk, depths, sample_names| {
        let depth_array = MultisampleDepthArray::from_samples(depths, sample_names.to_vec())?;
        let callable_mask =
            depth_array.apply_sample_thresholds(thresholds.min_depth, thresholds.max_depth);
        let callable_counts = depth_array.count_callable_per_population(
            &callable_mask,
            &pop_membership,
            thresholds.min_proportion,
            thresholds.mean_depth_range,
        )?;
        output_zarr.write_chunk(&chunk.contig_name, chunk.chunk_idx, callable_counts, None)?;
        Ok(())
    })?;

    Ok(output_zarr)
}

pub fn run_loci_zarr(
    input_zarr_path: PathBuf,
    output_path: PathBuf,
    pop_map: Option<PopulationMap>,
    thresholds: ThresholdConfig,
    output_per_sample_mask: bool,
) -> Result<()> {
    let spinner = create_spinner("Opening depth zarr...");
    let input_zarr = DepthArrays::open(&input_zarr_path)?;
    let input_zarr_samples = input_zarr.column_names();
    let chunk_size = input_zarr.chunk_size();

    spinner.set_message("Validating population map...");
    let pop_map =
        pop_map.unwrap_or_else(|| PopulationMap::default_from_samples(input_zarr_samples.to_vec()));
    pop_map.validate_exact_match(input_zarr_samples, false)?;
    spinner.finish_with_message(format!(
        "Loaded {} samples from zarr",
        input_zarr_samples.len()
    ));

    if output_per_sample_mask {
        process_sample_masks_zarr(input_zarr, output_path, thresholds, chunk_size, pop_map)?;
    } else {
        process_population_counts_zarr(input_zarr, output_path, pop_map, thresholds, chunk_size)?;
    }

    Ok(())
}

fn process_sample_masks_zarr(
    input_zarr: DepthArrays,
    output_path: PathBuf,
    thresholds: ThresholdConfig,
    chunk_size: u64,
    pop_map: PopulationMap,
) -> Result<SampleMaskArrays> {
    let output_zarr = SampleMaskArrays::create_new(
        output_path,
        input_zarr.contigs().clone(),
        input_zarr.column_names().to_vec(),
        chunk_size,
        Some(pop_map.populations_owned()),
    )?;

    let sample_names = input_zarr.column_names().to_vec();

    let total_chunks: usize = input_zarr
        .contigs()
        .iter()
        .map(|(_, len)| (len as u64).div_ceil(chunk_size) as usize)
        .sum();
    let pb = create_progress_bar(total_chunks);

    // Process by chromosome - open arrays once per chromosome
    for (chrom_name, chrom_length) in input_zarr.contigs().iter() {
        let num_chunks = (chrom_length as u64).div_ceil(chunk_size);

        let input_array = input_zarr.open_array(chrom_name)?;
        let output_array = output_zarr.open_array(chrom_name)?;

        (0..num_chunks).into_par_iter().try_for_each(|chunk_idx| {
            let depths = input_zarr.read_chunk(chrom_name, chunk_idx, Some(&input_array))?;

            let depth_array = MultisampleDepthArray {
                depths,
                sample_names: sample_names.clone(),
            };

            let callable_mask =
                depth_array.apply_sample_thresholds(thresholds.min_depth, thresholds.max_depth);

            output_zarr.write_chunk(chrom_name, chunk_idx, callable_mask, Some(&output_array))?;
            pb.inc(1);
            Ok::<_, color_eyre::Report>(())
        })?;
    }

    pb.finish_with_message("done");
    Ok(output_zarr)
}

fn process_population_counts_zarr(
    input_zarr: DepthArrays,
    output_path: PathBuf,
    pop_map: PopulationMap,
    thresholds: ThresholdConfig,
    chunk_size: u64,
) -> Result<CallableArrays> {
    let population_names: Vec<String> = pop_map.populations().map(|p| p.name.clone()).collect();

    let output_zarr = CallableArrays::create_new(
        output_path,
        input_zarr.contigs().clone(),
        population_names,
        chunk_size,
        Some(pop_map.populations_owned()),
    )?;

    let pop_membership = build_pop_membership(input_zarr.column_names(), &pop_map);
    let sample_names = input_zarr.column_names().to_vec();

    let total_chunks: usize = input_zarr
        .contigs()
        .iter()
        .map(|(_, len)| (len as u64).div_ceil(chunk_size) as usize)
        .sum();
    let pb = create_progress_bar(total_chunks);

    // Process by chromosome - open arrays once per chromosome
    for (chrom_name, chrom_length) in input_zarr.contigs().iter() {
        let num_chunks = (chrom_length as u64).div_ceil(chunk_size);

        let input_array = input_zarr.open_array(chrom_name)?;
        let output_array = output_zarr.open_array(chrom_name)?;

        (0..num_chunks).into_par_iter().try_for_each(|chunk_idx| {
            let depths = input_zarr.read_chunk(chrom_name, chunk_idx, Some(&input_array))?;

            let depth_array = MultisampleDepthArray {
                depths,
                sample_names: sample_names.clone(),
            };

            let callable_mask =
                depth_array.apply_sample_thresholds(thresholds.min_depth, thresholds.max_depth);

            let callable_counts = depth_array.count_callable_per_population(
                &callable_mask,
                &pop_membership,
                thresholds.min_proportion,
                thresholds.mean_depth_range,
            )?;

            output_zarr.write_chunk(chrom_name, chunk_idx, callable_counts, Some(&output_array))?;
            pb.inc(1);
            Ok::<_, color_eyre::Report>(())
        })?;
    }

    pb.finish_with_message("done");
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
        #[values((0.0, 1), (21.0, 0))] threshold_expected: (f64, u16),
    ) {
        let (min_depth, expected) = threshold_expected;
        let dir = TempDir::new().unwrap();
        let output_path = dir.path().join("callable.zarr");

        let depth_files = vec![PathBuf::from(path)];
        let config = SamplesConfig::from_paths(depth_files).unwrap();
        let pop_map = PopulationMap::default_from_samples(vec!["sample1".to_string()]);

        let thresholds = ThresholdConfig {
            min_depth,
            max_depth: f64::INFINITY,
            min_proportion: 0.0,
            mean_depth_range: (0.0, f64::INFINITY),
            per_contig: None,
        };

        run_loci(
            &config,
            output_path.clone(),
            Some(pop_map),
            thresholds,
            100,
            false,
            None,
        )
        .unwrap();

        assert!(output_path.exists());

        let arrays = CallableArrays::open(&output_path).unwrap();

        for (chrom_name, chrom_length) in arrays.contigs().iter() {
            let num_chunks = (chrom_length as u64 + 99) / 100;
            for chunk_idx in 0..num_chunks {
                let chunk_data = arrays.read_chunk(chrom_name, chunk_idx, None).unwrap();
                assert!(chunk_data.iter().all(|&x| x == expected));
            }
        }
    }

    #[rstest]
    fn test_two_samples(#[values((0.0, 2), (21.0, 0))] threshold_expected: (f64, u8)) {
        let (min_depth, _expected) = threshold_expected;
        let dir = TempDir::new().unwrap();
        let output_path = dir.path().join("callable.zarr");

        let depth_files = vec![
            PathBuf::from("tests/data/depth/all20/sample1.d4"),
            PathBuf::from("tests/data/depth/all20/sample2.d4"),
        ];
        let config = SamplesConfig::from_paths(depth_files).unwrap();

        let mut pop_data = indexmap::IndexMap::new();
        pop_data.insert("pop1".to_string(), vec!["sample1".to_string()]);
        pop_data.insert("pop2".to_string(), vec!["sample2".to_string()]);
        let pop_map = PopulationMap::from_populations(pop_data).unwrap();

        let thresholds = ThresholdConfig {
            min_depth,
            max_depth: f64::INFINITY,
            min_proportion: 0.0,
            mean_depth_range: (0.0, f64::INFINITY),
            per_contig: None,
        };

        run_loci(
            &config,
            output_path.clone(),
            Some(pop_map),
            thresholds,
            100,
            false,
            None,
        )
        .unwrap();

        assert!(output_path.exists());

        let arrays = CallableArrays::open(&output_path).unwrap();
        assert_eq!(arrays.column_names(), &["pop1", "pop2"]);

        for (chrom_name, chrom_length) in arrays.contigs().iter() {
            let num_chunks = (chrom_length as u64 + 99) / 100;
            for chunk_idx in 0..num_chunks {
                let chunk_data = arrays.read_chunk(chrom_name, chunk_idx, None).unwrap();

                let expected_rows = std::cmp::min(100, chrom_length - (chunk_idx * 100) as usize);
                assert_eq!(chunk_data.shape(), &[expected_rows, 2]);

                let expected_per_pop = if min_depth == 0.0 { 1 } else { 0 };
                assert!(chunk_data.iter().all(|&x| x == expected_per_pop));
            }
        }
    }

    #[rstest]
    fn test_per_sample_mask_output(
        #[values(
            "tests/data/depth/all20/sample1.g.vcf.gz",
            "tests/data/depth/all20/sample1.d4",
            "tests/data/depth/all20/sample1.d4.gz"
        )]
        path: &str,
        #[values((0.0, true), (21.0, false))] threshold_expected: (f64, bool),
    ) {
        let (min_depth, expected) = threshold_expected;
        let dir = TempDir::new().unwrap();
        let output_path = dir.path().join("masks.zarr");

        let depth_files = vec![PathBuf::from(path)];
        let config = SamplesConfig::from_paths(depth_files).unwrap();
        let pop_map = PopulationMap::default_from_samples(vec!["sample1".to_string()]);

        let thresholds = ThresholdConfig {
            min_depth,
            max_depth: f64::INFINITY,
            min_proportion: 0.0,
            mean_depth_range: (0.0, f64::INFINITY),
            per_contig: None,
        };

        run_loci(
            &config,
            output_path.clone(),
            Some(pop_map),
            thresholds,
            100,
            true,
            None,
        )
        .unwrap();

        assert!(output_path.exists());

        let arrays = SampleMaskArrays::open(&output_path).unwrap();
        assert_eq!(arrays.column_names(), &["sample1"]);

        for (chrom_name, chrom_length) in arrays.contigs().iter() {
            let num_chunks = (chrom_length as u64 + 99) / 100;
            for chunk_idx in 0..num_chunks {
                let chunk_data = arrays.read_chunk(chrom_name, chunk_idx, None).unwrap();
                assert!(chunk_data.iter().all(|&x| x == expected));
            }
        }
    }

    #[rstest]
    fn test_two_samples_per_sample_mask(
        #[values((0.0, true), (21.0, false))] threshold_expected: (f64, bool),
    ) {
        let (min_depth, expected) = threshold_expected;
        let dir = TempDir::new().unwrap();
        let output_path = dir.path().join("masks.zarr");

        let depth_files = vec![
            PathBuf::from("tests/data/depth/all20/sample1.d4"),
            PathBuf::from("tests/data/depth/all20/sample2.d4"),
        ];
        let config = SamplesConfig::from_paths(depth_files).unwrap();

        let mut pop_data = indexmap::IndexMap::new();
        pop_data.insert("pop1".to_string(), vec!["sample1".to_string()]);
        pop_data.insert("pop2".to_string(), vec!["sample2".to_string()]);
        let pop_map = PopulationMap::from_populations(pop_data).unwrap();

        let thresholds = ThresholdConfig {
            min_depth,
            max_depth: f64::INFINITY,
            min_proportion: 0.0,
            mean_depth_range: (0.0, f64::INFINITY),
            per_contig: None,
        };

        run_loci(
            &config,
            output_path.clone(),
            Some(pop_map),
            thresholds,
            100,
            true,
            None,
        )
        .unwrap();

        assert!(output_path.exists());

        let arrays = SampleMaskArrays::open(&output_path).unwrap();
        assert_eq!(arrays.column_names(), &["sample1", "sample2"]);

        for (chrom_name, chrom_length) in arrays.contigs().iter() {
            let num_chunks = (chrom_length as u64 + 99) / 100;
            for chunk_idx in 0..num_chunks {
                let chunk_data = arrays.read_chunk(chrom_name, chunk_idx, None).unwrap();

                let expected_rows = std::cmp::min(100, chrom_length - (chunk_idx * 100) as usize);
                assert_eq!(chunk_data.shape(), &[expected_rows, 2]);

                assert!(chunk_data.iter().all(|&x| x == expected));
            }
        }
    }

    /// A/B test: verify that SamplesConfig::from_file() produces same output as from_paths()
    #[test]
    fn test_samples_file_vs_from_paths_equivalence() {
        use std::io::Write;

        let dir = TempDir::new().unwrap();

        // Create samples TSV file
        let samples_tsv = dir.path().join("samples.tsv");
        {
            let mut f = std::fs::File::create(&samples_tsv).unwrap();
            writeln!(f, "sample_name\tfile_path\tpopulation").unwrap();
            writeln!(f, "sample1\ttests/data/depth/all20/sample1.d4\tpop1").unwrap();
            writeln!(f, "sample2\ttests/data/depth/all20/sample2.d4\tpop2").unwrap();
        }

        // Output paths
        let output_from_paths = dir.path().join("from_paths.zarr");
        let output_from_file = dir.path().join("from_file.zarr");

        let thresholds = ThresholdConfig {
            min_depth: 0.0,
            max_depth: f64::INFINITY,
            min_proportion: 0.0,
            mean_depth_range: (0.0, f64::INFINITY),
            per_contig: None,
        };

        // Method 1: from_paths with explicit pop_map
        let depth_files = vec![
            PathBuf::from("tests/data/depth/all20/sample1.d4"),
            PathBuf::from("tests/data/depth/all20/sample2.d4"),
        ];
        let config_from_paths = SamplesConfig::from_paths(depth_files).unwrap();

        let mut pop_data = indexmap::IndexMap::new();
        pop_data.insert("pop1".to_string(), vec!["sample1".to_string()]);
        pop_data.insert("pop2".to_string(), vec!["sample2".to_string()]);
        let pop_map = PopulationMap::from_populations(pop_data).unwrap();

        run_loci(
            &config_from_paths,
            output_from_paths.clone(),
            Some(pop_map),
            thresholds,
            100,
            false,
            None,
        )
        .unwrap();

        // Method 2: from_file (includes population in TSV)
        let config_from_file = SamplesConfig::from_file(&samples_tsv, true).unwrap();

        let thresholds2 = ThresholdConfig {
            min_depth: 0.0,
            max_depth: f64::INFINITY,
            min_proportion: 0.0,
            mean_depth_range: (0.0, f64::INFINITY),
            per_contig: None,
        };

        run_loci(
            &config_from_file,
            output_from_file.clone(),
            None, // Use population from config
            thresholds2,
            100,
            false,
            None,
        )
        .unwrap();

        // Compare outputs
        let result1 = CallableArrays::open(&output_from_paths).unwrap();
        let result2 = CallableArrays::open(&output_from_file).unwrap();

        assert_eq!(result1.column_names(), result2.column_names());
        // Compare contig lengths
        assert_eq!(result1.contigs().len(), result2.contigs().len());

        for (chrom_name, chrom_length) in result1.contigs().iter() {
            let num_chunks = (chrom_length as u64 + 99) / 100;
            for chunk_idx in 0..num_chunks {
                let data1 = result1.read_chunk(chrom_name, chunk_idx, None).unwrap();
                let data2 = result2.read_chunk(chrom_name, chunk_idx, None).unwrap();
                assert_eq!(data1, data2, "Mismatch at {}:{}", chrom_name, chunk_idx);
            }
        }
    }
}
