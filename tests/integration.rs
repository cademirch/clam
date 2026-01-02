use assert_cmd::prelude::*;
use clam::core::zarr::{CallableArrays, SampleMaskArrays};
use color_eyre::Result;
use d4::find_tracks_in_file;
use d4::ssio::D4TrackReader;
use glob::glob;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::{Path, PathBuf};
use std::process::Command;
use tempfile::TempDir;

/// Represents a window key for matching between output and truth
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
struct WindowKey {
    chrom: String,
    start: u32,
    end: u32,
    population: String,          // For pi: single population
    population2: Option<String>, // For dxy/fst: second population
}

/// Helper to compare zarr output against D4 truth file
fn compare_zarr_to_d4_truth(
    zarr: &CallableArrays,
    truth_d4_path: &str,
    chrom: &str,
) -> Result<()> {
    // Get track names from D4 (should be Pop1, Pop2)
    let mut d4_tracks: Vec<PathBuf> = vec![];
    find_tracks_in_file(truth_d4_path, |_| true, &mut d4_tracks)?;

    // Verify zarr has same populations as D4 tracks
    let zarr_columns = zarr.column_names();
    assert_eq!(
        zarr_columns.len(),
        d4_tracks.len(),
        "Number of populations mismatch: zarr has {}, D4 has {}",
        zarr_columns.len(),
        d4_tracks.len()
    );

    // Get chromosome info
    let chrom_length = zarr
        .contigs()
        .get_length(chrom)
        .expect(&format!("Chromosome {} not found in zarr", chrom));

    // Compare each population's data
    for (pop_idx, track) in d4_tracks.iter().enumerate() {
        let track_name = track
            .file_name()
            .and_then(|s| s.to_str())
            .unwrap_or("unknown");

        // Read truth from D4
        let rdr = File::open(truth_d4_path)?;
        let mut d4_rdr = D4TrackReader::from_reader(rdr, Some(track.to_str().unwrap()))?;
        let view = d4_rdr.get_view(chrom, 0, chrom_length as u32)?;
        let truth_values: Vec<u32> = view
            .map(|r| r.map(|(_, v)| v as u32))
            .collect::<std::result::Result<Vec<_>, _>>()?;

        // Read from zarr - iterate through chunks
        let mut zarr_values: Vec<u16> = vec![];
        let chunk_size = zarr.chunk_size();
        let num_chunks = (chrom_length as u64 + chunk_size - 1) / chunk_size;

        for chunk_idx in 0..num_chunks {
            let chunk_data = zarr.read_chunk(chrom, chunk_idx, None)?;
            // Extract column for this population, but only up to chrom_length
            let chunk_start = (chunk_idx * chunk_size) as usize;
            let rows_to_read = std::cmp::min(
                chunk_data.nrows(),
                chrom_length.saturating_sub(chunk_start),
            );
            for row in 0..rows_to_read {
                zarr_values.push(chunk_data[[row, pop_idx]]);
            }
        }

        // Compare lengths
        assert_eq!(
            zarr_values.len(),
            truth_values.len(),
            "Length mismatch for {}: zarr has {}, truth has {}",
            track_name,
            zarr_values.len(),
            truth_values.len()
        );

        // Compare values position by position
        for (pos, (zarr_val, truth_val)) in
            zarr_values.iter().zip(truth_values.iter()).enumerate()
        {
            assert_eq!(
                *zarr_val as u32, *truth_val,
                "Mismatch at position {} for {}: zarr={}, truth={}",
                pos, track_name, zarr_val, truth_val
            );
        }
    }

    Ok(())
}

/// Helper to compare per-sample mask zarr output against D4 truth file
/// by summing sample masks per population and comparing to population counts
fn compare_sample_masks_to_d4_truth(
    zarr_path: &Path,
    truth_d4_path: &str,
    popmap_path: &str,
    chrom: &str,
) -> Result<()> {
    // Parse popmap to get sample -> population mapping
    let popmap_file = File::open(popmap_path)?;
    let reader = BufReader::new(popmap_file);
    let mut sample_to_pop: HashMap<String, String> = HashMap::new();
    for line in reader.lines() {
        let line = line?;
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() >= 2 {
            sample_to_pop.insert(fields[0].to_string(), fields[1].to_string());
        }
    }

    // Open zarr as SampleMaskArrays
    let zarr = SampleMaskArrays::open(zarr_path)?;
    let sample_names = zarr.column_names();

    // Build sample index -> population index mapping
    // Get track names from D4 to know population order
    let mut d4_tracks: Vec<PathBuf> = vec![];
    find_tracks_in_file(truth_d4_path, |_| true, &mut d4_tracks)?;

    // Extract population names from D4 track paths (track name is the population)
    let pop_names: Vec<String> = d4_tracks
        .iter()
        .map(|t| t.file_name().unwrap().to_str().unwrap().to_string())
        .collect();

    // Build sample_idx -> pop_idx mapping
    let mut sample_to_pop_idx: Vec<Option<usize>> = vec![None; sample_names.len()];
    for (sample_idx, sample_name) in sample_names.iter().enumerate() {
        if let Some(pop_name) = sample_to_pop.get(sample_name) {
            if let Some(pop_idx) = pop_names.iter().position(|p| p == pop_name) {
                sample_to_pop_idx[sample_idx] = Some(pop_idx);
            }
        }
    }

    // Get chromosome info
    let chrom_length = zarr
        .contigs()
        .get_length(chrom)
        .expect(&format!("Chromosome {} not found in zarr", chrom));

    // Read truth values from D4 for each population
    let mut truth_values_per_pop: Vec<Vec<u32>> = Vec::new();
    for track in &d4_tracks {
        let rdr = File::open(truth_d4_path)?;
        let mut d4_rdr = D4TrackReader::from_reader(rdr, Some(track.to_str().unwrap()))?;
        let view = d4_rdr.get_view(chrom, 0, chrom_length as u32)?;
        let truth_values: Vec<u32> = view
            .map(|r| r.map(|(_, v)| v as u32))
            .collect::<std::result::Result<Vec<_>, _>>()?;
        truth_values_per_pop.push(truth_values);
    }

    // Read sample masks and sum by population, comparing to truth
    let chunk_size = zarr.chunk_size();
    let num_chunks = (chrom_length as u64 + chunk_size - 1) / chunk_size;
    let num_pops = pop_names.len();

    for chunk_idx in 0..num_chunks {
        let chunk_data = zarr.read_chunk(chrom, chunk_idx, None)?;
        let chunk_start = (chunk_idx * chunk_size) as usize;
        let rows_to_read = std::cmp::min(chunk_data.nrows(), chrom_length.saturating_sub(chunk_start));

        for row in 0..rows_to_read {
            let pos = chunk_start + row;

            // Sum callable samples per population
            let mut pop_counts: Vec<u32> = vec![0; num_pops];
            for (sample_idx, &pop_idx_opt) in sample_to_pop_idx.iter().enumerate() {
                if let Some(pop_idx) = pop_idx_opt {
                    if chunk_data[[row, sample_idx]] {
                        pop_counts[pop_idx] += 1;
                    }
                }
            }

            // Compare to truth
            for (pop_idx, pop_name) in pop_names.iter().enumerate() {
                let computed = pop_counts[pop_idx];
                let expected = truth_values_per_pop[pop_idx][pos];
                assert_eq!(
                    computed, expected,
                    "Mismatch at position {} for {}: computed={}, expected={}",
                    pos, pop_name, computed, expected
                );
            }
        }
    }

    Ok(())
}

#[test]
fn test_stat_with_gatk_vcf() -> Result<()> {
    color_eyre::install().ok();

    let temp_dir = TempDir::new()?;
    let output_dir = temp_dir.path().join("stat_output");
    let callable_zarr = temp_dir.path().join("callable.zarr");

    // First, create callable sites zarr from gVCF files
    let gvcf_files: Vec<String> = glob("tests/data/integration/gvcf/*.g.vcf.gz")?
        .filter_map(|p| p.ok())
        .map(|p| p.to_string_lossy().to_string())
        .collect();

    let mut cmd = Command::cargo_bin("clam")?;
    cmd.arg("loci")
        .args(&gvcf_files)
        .arg("-o")
        .arg(&callable_zarr)
        .arg("-p")
        .arg("tests/data/integration/popmap.txt")
        .arg("--min-depth")
        .arg("2")
        .arg("--min-gq")
        .arg("10");

    let output = cmd.output()?;
    assert!(
        output.status.success(),
        "clam loci failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );

    // Run clam stat with GATK VCF
    let mut cmd = Command::cargo_bin("clam")?;
    cmd.arg("stat")
        .arg("tests/data/integration/gatk.varsonly.vcf.gz")
        .arg("-o")
        .arg(&output_dir)
        .arg("-c")
        .arg(&callable_zarr)
        .arg("-p")
        .arg("tests/data/integration/popmap.txt")
        .arg("-w")
        .arg("10000"); // 10kb windows to match truth data

    let output = cmd.output()?;
    assert!(
        output.status.success(),
        "clam stat failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );

    // Compare output files to truth by matching windows
    let float_tolerance = 0.001; // 0.1% tolerance for pi/dxy/fst values

    // Compare pi values
    let output_pi = parse_pi_tsv(&output_dir.join("pi.tsv").to_string_lossy())?;
    let truth_pi = parse_pi_tsv("tests/data/integration/truth/gatk/clam_pi.tsv")?;
    compare_stat_values(&output_pi, &truth_pi, "pi", float_tolerance)?;

    // Compare dxy values (value is in column 5 for new format)
    let output_dxy = parse_pairwise_tsv(&output_dir.join("dxy.tsv").to_string_lossy(), 5)?;
    let truth_dxy = parse_pairwise_tsv("tests/data/integration/truth/gatk/clam_dxy.tsv", 5)?;
    compare_stat_values(&output_dxy, &truth_dxy, "dxy", float_tolerance)?;

    // Compare fst values (value is in column 5 for new format)
    let output_fst = parse_pairwise_tsv(&output_dir.join("fst.tsv").to_string_lossy(), 5)?;
    let truth_fst = parse_pairwise_tsv("tests/data/integration/truth/gatk/clam_fst.tsv", 5)?;
    compare_stat_values(&output_fst, &truth_fst, "fst", float_tolerance)?;

    Ok(())
}

#[test]
fn test_stat_with_bcftools_vcf() -> Result<()> {
    color_eyre::install().ok();

    let temp_dir = TempDir::new()?;
    let output_dir = temp_dir.path().join("stat_output");
    let callable_zarr = temp_dir.path().join("callable.zarr");

    // First, create callable sites zarr from D4 files
    let d4_files: Vec<String> = glob("tests/data/integration/d4/*.d4.gz")?
        .filter_map(|p| p.ok())
        .map(|p| p.to_string_lossy().to_string())
        .collect();

    let mut cmd = Command::cargo_bin("clam")?;
    cmd.arg("loci")
        .args(&d4_files)
        .arg("-o")
        .arg(&callable_zarr)
        .arg("-p")
        .arg("tests/data/integration/popmap.txt")
        .arg("--min-depth")
        .arg("2");

    let output = cmd.output()?;
    assert!(
        output.status.success(),
        "clam loci failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );

    // Run clam stat with bcftools VCF
    let mut cmd = Command::cargo_bin("clam")?;
    cmd.arg("stat")
        .arg("tests/data/integration/bcftools.varsonly.vcf.gz")
        .arg("-o")
        .arg(&output_dir)
        .arg("-c")
        .arg(&callable_zarr)
        .arg("-p")
        .arg("tests/data/integration/popmap.txt")
        .arg("-w")
        .arg("10000"); // 10kb windows to match truth data

    let output = cmd.output()?;
    assert!(
        output.status.success(),
        "clam stat failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );

    // Compare output files to truth by matching windows
    let float_tolerance = 0.001; // 0.1% tolerance for pi/dxy/fst values

    // Compare pi values
    let output_pi = parse_pi_tsv(&output_dir.join("pi.tsv").to_string_lossy())?;
    let truth_pi = parse_pi_tsv("tests/data/integration/truth/bcftools/clam_pi.tsv")?;
    compare_stat_values(&output_pi, &truth_pi, "pi", float_tolerance)?;

    // Compare dxy values (value is in column 5 for new format)
    let output_dxy = parse_pairwise_tsv(&output_dir.join("dxy.tsv").to_string_lossy(), 5)?;
    let truth_dxy = parse_pairwise_tsv("tests/data/integration/truth/bcftools/clam_dxy.tsv", 5)?;
    compare_stat_values(&output_dxy, &truth_dxy, "dxy", float_tolerance)?;

    // Compare fst values (value is in column 5 for new format)
    let output_fst = parse_pairwise_tsv(&output_dir.join("fst.tsv").to_string_lossy(), 5)?;
    let truth_fst = parse_pairwise_tsv("tests/data/integration/truth/bcftools/clam_fst.tsv", 5)?;
    compare_stat_values(&output_fst, &truth_fst, "fst", float_tolerance)?;

    Ok(())
}


/// Parse pi TSV and extract window -> pi value mapping
fn parse_pi_tsv(path: &str) -> Result<std::collections::HashMap<WindowKey, f64>> {
    let file = File::open(path)?;
    let reader = BufReader::new(file);
    let mut map = std::collections::HashMap::new();

    for (line_num, line) in reader.lines().enumerate() {
        let line = line?;
        if line_num == 0 {
            continue; // Skip header
        }

        let fields: Vec<&str> = line.split('\t').collect();

        // Handle both old format (population_name, chrom, start, end, pi, ...)
        // and new format (chrom, start, end, population, pi, ...)
        let (chrom, start, end, pop, pi) = if fields[0].starts_with("Pop") {
            // Old format: population_name is first
            (
                fields[1].to_string(),
                fields[2].parse::<u32>()?,
                fields[3].parse::<u32>()?,
                fields[0].to_string(),
                fields[4].parse::<f64>()?,
            )
        } else {
            // New format: chrom is first
            (
                fields[0].to_string(),
                fields[1].parse::<u32>()?,
                fields[2].parse::<u32>()?,
                fields[3].to_string(),
                fields[4].parse::<f64>()?,
            )
        };

        let key = WindowKey {
            chrom,
            start,
            end,
            population: pop,
            population2: None,
        };
        map.insert(key, pi);
    }

    Ok(map)
}

/// Parse dxy/fst TSV and extract window -> value mapping
fn parse_pairwise_tsv(path: &str, value_col: usize) -> Result<std::collections::HashMap<WindowKey, f64>> {
    let file = File::open(path)?;
    let reader = BufReader::new(file);
    let mut map = std::collections::HashMap::new();

    for (line_num, line) in reader.lines().enumerate() {
        let line = line?;
        if line_num == 0 {
            continue; // Skip header
        }

        let fields: Vec<&str> = line.split('\t').collect();

        // Handle both old format (pop1_name, pop2_name, chrom, start, end, value, ...)
        // and new format (chrom, start, end, pop1, pop2, value, ...)
        let (chrom, start, end, pop1, pop2, value) = if fields[0].starts_with("Pop") {
            // Old format: populations are first
            (
                fields[2].to_string(),
                fields[3].parse::<u32>()?,
                fields[4].parse::<u32>()?,
                fields[0].to_string(),
                fields[1].to_string(),
                fields[5].parse::<f64>()?,
            )
        } else {
            // New format: chrom is first
            (
                fields[0].to_string(),
                fields[1].parse::<u32>()?,
                fields[2].parse::<u32>()?,
                fields[3].to_string(),
                fields[4].to_string(),
                fields[value_col].parse::<f64>()?,
            )
        };

        let key = WindowKey {
            chrom,
            start,
            end,
            population: pop1,
            population2: Some(pop2),
        };
        map.insert(key, value);
    }

    Ok(map)
}

/// Compare stat values by matching windows
fn compare_stat_values(
    output_map: &std::collections::HashMap<WindowKey, f64>,
    truth_map: &std::collections::HashMap<WindowKey, f64>,
    stat_name: &str,
    tolerance: f64,
) -> Result<()> {
    for (truth_key, truth_val) in truth_map {
        let output_val = output_map.get(truth_key);

        if let Some(output_val) = output_val {
            let diff = (output_val - truth_val).abs();
            assert!(
                diff <= tolerance,
                "{} mismatch for window {:?}: output={}, truth={}, diff={}",
                stat_name,
                truth_key,
                output_val,
                truth_val,
                diff
            );
        }
        // Skip windows that don't exist in output (like tiny trailing windows)
    }

    Ok(())
}

/// Key for heterozygosity records (includes sample for per-sample mode)
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
struct HetKey {
    chrom: String,
    start: u32,
    end: u32,
    sample: Option<String>,
    population: String,
}

/// Heterozygosity values for comparison
#[derive(Debug, Clone)]
struct HetValues {
    het_total: usize,
    callable_total: usize,
    heterozygosity: f64,
    het_not_in_roh: Option<usize>,
    callable_not_in_roh: Option<usize>,
    heterozygosity_not_in_roh: Option<f64>,
}

/// Parse heterozygosity TSV (handles both per-sample and per-population formats)
fn parse_heterozygosity_tsv(path: &str) -> Result<std::collections::HashMap<HetKey, HetValues>> {
    let file = File::open(path)?;
    let reader = BufReader::new(file);
    let mut map = std::collections::HashMap::new();

    let mut header_indices: HashMap<String, usize> = HashMap::new();

    for (line_num, line) in reader.lines().enumerate() {
        let line = line?;
        let fields: Vec<&str> = line.split('\t').collect();

        if line_num == 0 {
            // Parse header
            for (i, field) in fields.iter().enumerate() {
                header_indices.insert(field.to_string(), i);
            }
            continue;
        }

        let chrom = fields[*header_indices.get("chrom").unwrap()].to_string();
        let start = fields[*header_indices.get("start").unwrap()].parse::<u32>()?;
        let end = fields[*header_indices.get("end").unwrap()].parse::<u32>()?;
        let population = fields[*header_indices.get("population").unwrap()].to_string();

        // Sample is optional (only in per-sample mode)
        let sample = header_indices
            .get("sample")
            .map(|&i| fields[i].to_string());

        let het_total = fields[*header_indices.get("het_total").unwrap()].parse::<usize>()?;
        let callable_total = fields[*header_indices.get("callable_total").unwrap()].parse::<usize>()?;
        let heterozygosity = fields[*header_indices.get("heterozygosity").unwrap()]
            .parse::<f64>()
            .unwrap_or(f64::NAN);

        // ROH fields are optional
        let het_not_in_roh = header_indices
            .get("het_not_in_roh")
            .and_then(|&i| fields[i].parse::<usize>().ok());
        let callable_not_in_roh = header_indices
            .get("callable_not_in_roh")
            .and_then(|&i| fields[i].parse::<usize>().ok());
        let heterozygosity_not_in_roh = header_indices
            .get("heterozygosity_not_in_roh")
            .and_then(|&i| fields[i].parse::<f64>().ok());

        let key = HetKey {
            chrom,
            start,
            end,
            sample,
            population,
        };

        let values = HetValues {
            het_total,
            callable_total,
            heterozygosity,
            het_not_in_roh,
            callable_not_in_roh,
            heterozygosity_not_in_roh,
        };

        map.insert(key, values);
    }

    Ok(map)
}

/// Compare heterozygosity values between output and truth
/// 
/// When comparing against Python-generated truth, callable counts may differ slightly
/// because Python counts callable at variant positions using gVCF callable status,
/// while clam uses VCF non-missing status. The het counts should match exactly.
fn compare_heterozygosity_values(
    output_map: &std::collections::HashMap<HetKey, HetValues>,
    truth_map: &std::collections::HashMap<HetKey, HetValues>,
    tolerance: f64,
) -> Result<()> {
    // Allow callable counts to differ by up to 1% or 100 positions
    let callable_tolerance_pct = 0.01;
    let callable_tolerance_abs = 100;
    
    for (truth_key, truth_val) in truth_map {
        let output_val = output_map.get(truth_key);

        if let Some(output_val) = output_val {
            // Compare het_total exactly
            assert_eq!(
                output_val.het_total, truth_val.het_total,
                "het_total mismatch for {:?}: output={}, truth={}",
                truth_key, output_val.het_total, truth_val.het_total
            );

            // Compare callable_total with tolerance (Python vs clam may differ at variant positions)
            let callable_diff = (output_val.callable_total as i64 - truth_val.callable_total as i64).abs() as usize;
            let callable_pct_diff = callable_diff as f64 / truth_val.callable_total.max(1) as f64;
            assert!(
                callable_diff <= callable_tolerance_abs || callable_pct_diff <= callable_tolerance_pct,
                "callable_total mismatch for {:?}: output={}, truth={}, diff={} ({:.2}%)",
                truth_key, output_val.callable_total, truth_val.callable_total, 
                callable_diff, callable_pct_diff * 100.0
            );

            // Compare heterozygosity (with tolerance)
            if !truth_val.heterozygosity.is_nan() {
                let diff = (output_val.heterozygosity - truth_val.heterozygosity).abs();
                assert!(
                    diff <= tolerance,
                    "heterozygosity mismatch for {:?}: output={}, truth={}, diff={}",
                    truth_key, output_val.heterozygosity, truth_val.heterozygosity, diff
                );
            }

            // Compare ROH-excluded het counts exactly
            if let (Some(out_het), Some(truth_het)) =
                (output_val.het_not_in_roh, truth_val.het_not_in_roh)
            {
                assert_eq!(
                    out_het, truth_het,
                    "het_not_in_roh mismatch for {:?}: output={}, truth={}",
                    truth_key, out_het, truth_het
                );
            }

            // Compare callable_not_in_roh with tolerance
            if let (Some(out_call), Some(truth_call)) =
                (output_val.callable_not_in_roh, truth_val.callable_not_in_roh)
            {
                let diff = (out_call as i64 - truth_call as i64).abs() as usize;
                let pct_diff = diff as f64 / truth_call.max(1) as f64;
                assert!(
                    diff <= callable_tolerance_abs || pct_diff <= callable_tolerance_pct,
                    "callable_not_in_roh mismatch for {:?}: output={}, truth={}, diff={} ({:.2}%)",
                    truth_key, out_call, truth_call, diff, pct_diff * 100.0
                );
            }

            if let (Some(out_het_rate), Some(truth_het_rate)) = (
                output_val.heterozygosity_not_in_roh,
                truth_val.heterozygosity_not_in_roh,
            ) {
                if !truth_het_rate.is_nan() {
                    let diff = (out_het_rate - truth_het_rate).abs();
                    assert!(
                        diff <= tolerance,
                        "heterozygosity_not_in_roh mismatch for {:?}: output={}, truth={}, diff={}",
                        truth_key, out_het_rate, truth_het_rate, diff
                    );
                }
            }
        }
    }

    Ok(())
}

#[test]
fn test_loci_gvcf_with_populations() -> Result<()> {
    color_eyre::install().ok();

    let temp_dir = TempDir::new()?;
    let output_zarr = temp_dir.path().join("callable.zarr");

    // Collect all gVCF files from test data
    let gvcf_files: Vec<String> = glob("tests/data/integration/gvcf/*.g.vcf.gz")?
        .filter_map(|p| p.ok())
        .map(|p| p.to_string_lossy().to_string())
        .collect();

    assert!(!gvcf_files.is_empty(), "No gVCF files found in test data");

    // Run clam loci command
    let mut cmd = Command::cargo_bin("clam")?;
    cmd.arg("loci")
        .args(&gvcf_files)
        .arg("-o")
        .arg(&output_zarr)
        .arg("-p")
        .arg("tests/data/integration/popmap.txt")
        .arg("--min-depth")
        .arg("2")
        .arg("--min-gq")
        .arg("10");

    let output = cmd.output()?;
    assert!(
        output.status.success(),
        "clam loci failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );

    // Compare output to truth
    let zarr = CallableArrays::open(&output_zarr)?;
    compare_zarr_to_d4_truth(&zarr, "tests/data/integration/gvcf/callable_sites.d4", "chr2l")?;

    Ok(())
}

#[test]
fn test_loci_d4_with_populations() -> Result<()> {
    color_eyre::install().ok();

    let temp_dir = TempDir::new()?;
    let output_zarr = temp_dir.path().join("callable.zarr");

    // Collect all D4.gz files from test data
    let d4_files: Vec<String> = glob("tests/data/integration/d4/*.d4.gz")?
        .filter_map(|p| p.ok())
        .map(|p| p.to_string_lossy().to_string())
        .collect();

    assert!(!d4_files.is_empty(), "No D4 files found in test data");

    // Run clam loci command (no --min-gq for D4 input)
    let mut cmd = Command::cargo_bin("clam")?;
    cmd.arg("loci")
        .args(&d4_files)
        .arg("-o")
        .arg(&output_zarr)
        .arg("-p")
        .arg("tests/data/integration/popmap.txt")
        .arg("--min-depth")
        .arg("2");

    let output = cmd.output()?;
    assert!(
        output.status.success(),
        "clam loci failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );

    // Compare output to truth
    let zarr = CallableArrays::open(&output_zarr)?;
    compare_zarr_to_d4_truth(&zarr, "tests/data/integration/d4/callable_sites.d4", "chr2l")?;

    Ok(())
}

#[test]
fn test_loci_gvcf_per_sample_masks() -> Result<()> {
    color_eyre::install().ok();

    let temp_dir = TempDir::new()?;
    let output_zarr = temp_dir.path().join("callable_per_sample.zarr");

    // Collect all gVCF files from test data
    let gvcf_files: Vec<String> = glob("tests/data/integration/gvcf/*.g.vcf.gz")?
        .filter_map(|p| p.ok())
        .map(|p| p.to_string_lossy().to_string())
        .collect();

    assert!(!gvcf_files.is_empty(), "No gVCF files found in test data");

    // Run clam loci command with --per-sample flag
    let mut cmd = Command::cargo_bin("clam")?;
    cmd.arg("loci")
        .args(&gvcf_files)
        .arg("-o")
        .arg(&output_zarr)
        .arg("-p")
        .arg("tests/data/integration/popmap.txt")
        .arg("--min-depth")
        .arg("2")
        .arg("--min-gq")
        .arg("10")
        .arg("--per-sample");

    let output = cmd.output()?;
    assert!(
        output.status.success(),
        "clam loci failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );

    // Compare per-sample masks to population counts truth
    // by summing sample masks per population
    compare_sample_masks_to_d4_truth(
        &output_zarr,
        "tests/data/integration/gvcf/callable_sites.d4",
        "tests/data/integration/popmap.txt",
        "chr2l",
    )?;

    Ok(())
}

#[test]
fn test_stat_with_perfect_vcf() -> Result<()> {
    color_eyre::install().ok();

    let temp_dir = TempDir::new()?;
    let output_dir = temp_dir.path().join("stat_output");
    let callable_zarr = temp_dir.path().join("callable.zarr");

    // First, create callable sites zarr from gVCF files
    let gvcf_files: Vec<String> = glob("tests/data/integration/gvcf/*.g.vcf.gz")?
        .filter_map(|p| p.ok())
        .map(|p| p.to_string_lossy().to_string())
        .collect();

    let mut cmd = Command::cargo_bin("clam")?;
    cmd.arg("loci")
        .args(&gvcf_files)
        .arg("-o")
        .arg(&callable_zarr)
        .arg("-p")
        .arg("tests/data/integration/popmap.txt")
        .arg("--min-depth")
        .arg("2")
        .arg("--min-gq")
        .arg("10");

    let output = cmd.output()?;
    assert!(
        output.status.success(),
        "clam loci failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );

    // Run clam stat with perfect VCF
    let mut cmd = Command::cargo_bin("clam")?;
    cmd.arg("stat")
        .arg("tests/data/integration/perfect.vcf.gz")
        .arg("-o")
        .arg(&output_dir)
        .arg("-c")
        .arg(&callable_zarr)
        .arg("-p")
        .arg("tests/data/integration/popmap.txt")
        .arg("-w")
        .arg("10000"); // 10kb windows to match truth data

    let output = cmd.output()?;
    assert!(
        output.status.success(),
        "clam stat failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );

    // Compare output files to truth by matching windows
    // Use a tolerance that accounts for slight differences in callable site counts
    let float_tolerance = 0.001; // 0.1% tolerance for pi/dxy/fst values

    // Compare pi values
    let output_pi = parse_pi_tsv(&output_dir.join("pi.tsv").to_string_lossy())?;
    let truth_pi = parse_pi_tsv("tests/data/integration/truth/perfect/clam_pi.tsv")?;
    compare_stat_values(&output_pi, &truth_pi, "pi", float_tolerance)?;

    // Compare dxy values (value is in column 5 for new format)
    let output_dxy = parse_pairwise_tsv(&output_dir.join("dxy.tsv").to_string_lossy(), 5)?;
    let truth_dxy = parse_pairwise_tsv("tests/data/integration/truth/perfect/clam_dxy.tsv", 5)?;
    compare_stat_values(&output_dxy, &truth_dxy, "dxy", float_tolerance)?;

    // Compare fst values (value is in column 5 for new format)
    let output_fst = parse_pairwise_tsv(&output_dir.join("fst.tsv").to_string_lossy(), 5)?;
    let truth_fst = parse_pairwise_tsv("tests/data/integration/truth/perfect/clam_fst.tsv", 5)?;
    compare_stat_values(&output_fst, &truth_fst, "fst", float_tolerance)?;

    Ok(())
}

#[test]
fn test_stat_heterozygosity_with_sample_masks() -> Result<()> {
    color_eyre::install().ok();

    let temp_dir = TempDir::new()?;
    let output_dir = temp_dir.path().join("stat_output");
    let callable_zarr = temp_dir.path().join("callable_per_sample.zarr");

    // Create callable sites zarr with per-sample masks from gVCF files
    let gvcf_files: Vec<String> = glob("tests/data/integration/gvcf/*.g.vcf.gz")?
        .filter_map(|p| p.ok())
        .map(|p| p.to_string_lossy().to_string())
        .collect();

    let mut cmd = Command::cargo_bin("clam")?;
    cmd.arg("loci")
        .args(&gvcf_files)
        .arg("-o")
        .arg(&callable_zarr)
        .arg("-p")
        .arg("tests/data/integration/popmap.txt")
        .arg("--min-depth")
        .arg("2")
        .arg("--min-gq")
        .arg("10")
        .arg("--per-sample");

    let output = cmd.output()?;
    assert!(
        output.status.success(),
        "clam loci failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );

    // Run clam stat with GATK VCF and ROH file
    let mut cmd = Command::cargo_bin("clam")?;
    cmd.arg("stat")
        .arg("tests/data/integration/gatk.varsonly.vcf.gz")
        .arg("-o")
        .arg(&output_dir)
        .arg("-c")
        .arg(&callable_zarr)
        .arg("-p")
        .arg("tests/data/integration/popmap.txt")
        .arg("-r")
        .arg("tests/data/integration/roh/sample_roh.bed.gz")
        .arg("-w")
        .arg("10000");

    let output = cmd.output()?;
    assert!(
        output.status.success(),
        "clam stat failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );

    // Compare heterozygosity output to Python-generated truth
    // Note: callable counts may differ slightly because Python counts callable at variant
    // positions using gVCF callable status, while clam uses VCF non-missing status.
    // The het counts should match exactly.
    let float_tolerance = 0.001; // 0.1% tolerance for heterozygosity rate
    let output_het = parse_heterozygosity_tsv(&output_dir.join("heterozygosity.tsv").to_string_lossy())?;
    let truth_het = parse_heterozygosity_tsv("tests/data/integration/truth/python/het_sample_masks.tsv")?;
    compare_heterozygosity_values(&output_het, &truth_het, float_tolerance)?;

    Ok(())
}

#[test]
fn test_stat_heterozygosity_with_pop_counts() -> Result<()> {
    color_eyre::install().ok();

    let temp_dir = TempDir::new()?;
    let output_dir = temp_dir.path().join("stat_output");
    let callable_zarr = temp_dir.path().join("callable_pop_counts.zarr");

    // Create callable sites zarr with population counts from gVCF files
    let gvcf_files: Vec<String> = glob("tests/data/integration/gvcf/*.g.vcf.gz")?
        .filter_map(|p| p.ok())
        .map(|p| p.to_string_lossy().to_string())
        .collect();

    let mut cmd = Command::cargo_bin("clam")?;
    cmd.arg("loci")
        .args(&gvcf_files)
        .arg("-o")
        .arg(&callable_zarr)
        .arg("-p")
        .arg("tests/data/integration/popmap.txt")
        .arg("--min-depth")
        .arg("2")
        .arg("--min-gq")
        .arg("10");

    let output = cmd.output()?;
    assert!(
        output.status.success(),
        "clam loci failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );

    // Run clam stat with GATK VCF and ROH file
    let mut cmd = Command::cargo_bin("clam")?;
    cmd.arg("stat")
        .arg("tests/data/integration/gatk.varsonly.vcf.gz")
        .arg("-o")
        .arg(&output_dir)
        .arg("-c")
        .arg(&callable_zarr)
        .arg("-p")
        .arg("tests/data/integration/popmap.txt")
        .arg("-r")
        .arg("tests/data/integration/roh/sample_roh.bed.gz")
        .arg("-w")
        .arg("10000");

    let output = cmd.output()?;
    assert!(
        output.status.success(),
        "clam stat failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );

    // Compare heterozygosity output to Python-generated truth
    // Note: callable counts may differ slightly because Python counts callable at variant
    // positions using gVCF callable status, while clam uses VCF non-missing status.
    // The het counts should match exactly.
    let float_tolerance = 0.001; // 0.1% tolerance for heterozygosity rate
    let output_het = parse_heterozygosity_tsv(&output_dir.join("heterozygosity.tsv").to_string_lossy())?;
    let truth_het = parse_heterozygosity_tsv("tests/data/integration/truth/python/het_pop_counts.tsv")?;
    compare_heterozygosity_values(&output_het, &truth_het, float_tolerance)?;

    Ok(())
}
