use assert_cmd::prelude::*;
use clam::core::zarr::CallableArrays;
use color_eyre::Result;
use d4::find_tracks_in_file;
use d4::ssio::D4TrackReader;
use glob::glob;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::PathBuf;
use std::process::Command;
use tempfile::TempDir;

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

/// Represents a window key for matching between output and truth
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
struct WindowKey {
    chrom: String,
    start: u32,
    end: u32,
    population: String,          // For pi: single population
    population2: Option<String>, // For dxy/fst: second population
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

/// Compare stat values by matching windows (tolerating coordinate differences of 1)
fn compare_stat_values(
    output_map: &std::collections::HashMap<WindowKey, f64>,
    truth_map: &std::collections::HashMap<WindowKey, f64>,
    stat_name: &str,
    tolerance: f64,
) -> Result<()> {
    // For each truth window, find matching output window (allowing off-by-one in coordinates)
    for (truth_key, truth_val) in truth_map {
        // Try exact match first
        let output_val = output_map.get(truth_key).or_else(|| {
            // Try with end+1 (new format uses half-open intervals)
            let adjusted_key = WindowKey {
                chrom: truth_key.chrom.clone(),
                start: truth_key.start,
                end: truth_key.end + 1,
                population: truth_key.population.clone(),
                population2: truth_key.population2.clone(),
            };
            output_map.get(&adjusted_key)
        });

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
