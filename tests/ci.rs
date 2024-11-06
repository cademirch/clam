use anyhow::{Context, Result};
use assert_cmd::prelude::*;
use log::debug;
use std::fs;
use std::io::{self, BufReader, Read, BufRead};
use std::path::{Path, PathBuf};
use std::process::Command;
use std::str::FromStr;
use std::fs::File;
use d4::ssio::{D4MatrixReader, D4TrackReader};
use d4::find_tracks_in_file;

/// Compare two files byte-by-byte to check if they are the same.
fn files_are_equal<P: AsRef<Path> + Clone>(file1: P, file2: P) -> Result<bool> {
    let f1 = File::open(file1.clone()).context(format!("{} does not exist", file1.as_ref().display()))?;
    let f2 = File::open(file2.clone()).context(format!("{} does not exist", file2.as_ref().display()))?;

    let mut reader1 = BufReader::new(f1);
    let mut reader2 = BufReader::new(f2);

    let mut buffer1 = [0; 1024];
    let mut buffer2 = [0; 1024];

    loop {
        let bytes_read1 = reader1.read(&mut buffer1)?;
        let bytes_read2 = reader2.read(&mut buffer2)?;

        if bytes_read1 != bytes_read2 {
            return Ok(false); // Files are different sizes
        }

        if bytes_read1 == 0 {
            break; // End of file reached and files are the same
        }

        if buffer1[..bytes_read1] != buffer2[..bytes_read2] {
            return Ok(false); // Content mismatch
        }
    }

    Ok(true)
}

fn d4_files_are_equal<P: AsRef<Path> + Clone>(file1: P, file2: P, chrom: &str, begin: u32, end: u32) -> Result<bool> {
    let mut f1_tracks: Vec<PathBuf> = vec![];
    find_tracks_in_file(file1.clone(), |_| true, &mut f1_tracks)?;

    let mut f2_tracks: Vec<PathBuf> = vec![];
    find_tracks_in_file(file2.clone(), |_| true, &mut f2_tracks)?;

    if f1_tracks.len() != f2_tracks.len() {
        return Ok(false); // Different number of tracks means files are not equal
    }

    let mut f1_data: Vec<Vec<u32>> = vec![vec![]; f1_tracks.len()];
    let mut f2_data: Vec<Vec<u32>> = vec![vec![]; f2_tracks.len()];

    // Read data for each track in f1
    for (idx, track) in f1_tracks.iter().enumerate() {
        let rdr = File::open(file1.clone())?;
        let mut d4_rdr = D4TrackReader::from_reader(rdr, Some(track.to_str().unwrap()))?;
        let view = d4_rdr.get_view(chrom, begin, end)?;

        for result in view {
            let (pos, value) = result?;
            f1_data[idx].push(value as u32);  // Assuming values are within u32 range
        }
    }

    // Read data for each track in f2
    for (idx, track) in f2_tracks.iter().enumerate() {
        let rdr = File::open(file2.clone())?;
        let mut d4_rdr = D4TrackReader::from_reader(rdr, Some(track.to_str().unwrap()))?;
        let view = d4_rdr.get_view(chrom, begin, end)?;

        for result in view {
            let (pos, value) = result?;
            f2_data[idx].push(value as u32);
        }
    }

    // Compare data from f1 and f2
    for (track1, track2) in f1_data.iter().zip(f2_data.iter()) {
        if track1 != track2 {
            return Ok(false);  // Mismatch found
        }
    }

    Ok(true)  // All data matched
}




#[test]
fn loci_test_no_pops_output_bed() -> Result<(), Box<dyn std::error::Error>> {
    //bamsim generate --num-bams 5 --num-chrs 1 --chr-length 1000 --min-depth 3 --max-depth 15 --min-mean-depth 6 --proportion 0.5 test_no_pops/
    let _ = env_logger::builder()
        .target(env_logger::Target::Stdout)
        .filter_level(log::LevelFilter::Trace)
        .is_test(true)
        .try_init();
    let mut cmd = Command::cargo_bin("clam")?;

    let input_file = PathBuf::from_str("tests/data/test_no_pops/merged.d4")?;
    let output_prefix = PathBuf::from_str("tests/data/test_no_pops/output_bed")?;
    let output_file = output_prefix.with_extension("bed");
    let expected_output_file = PathBuf::from_str("tests/data/test_no_pops/truth.bed")?;

    cmd.arg("loci")
        .arg(&input_file) // infile as a positional argument
        .arg(&output_prefix) // unique outfile as a positional argument
        .arg("-m")
        .arg("3")
        .arg("-M")
        .arg("15")
        .arg("-d")
        .arg("0.5")
        .arg("-u")
        .arg("6")
        .arg("--no-counts");

    let output = cmd.output()?;
    debug!("Ran command");
    debug!("stdout: {}", String::from_utf8_lossy(&output.stdout));
    debug!("stderr: {}", String::from_utf8_lossy(&output.stderr));

    let success = files_are_equal(expected_output_file.clone(), output_file.clone())?;
    assert_eq!(true, success);
    fs::remove_file(&output_file)?;

    Ok(())
}

#[test]
fn loci_test_pops() -> Result<(), Box<dyn std::error::Error>> {
    //./target/release/bamsim generate --num-bams 5 --num-chrs 1 --chr-length 1000 --min-depth 3 --max-depth 15 --min-mean-depth 6 --proportion 0.5 test_pops/pop1/bams/
    let _ = env_logger::builder()
        .target(env_logger::Target::Stdout)
        .filter_level(log::LevelFilter::Trace)
        .is_test(true)
        .try_init();
    let mut cmd = Command::cargo_bin("clam")?;

    let input_file = PathBuf::from_str("tests/data/test_pops/merged.d4")?;
    let output_prefix = PathBuf::from_str("tests/data/test_pops/output_d4")?;
    let output_file = output_prefix.with_extension("d4");
    let expected_output_file = PathBuf::from_str("tests/data/test_pops/counts.d4")?;
    let pop_file = PathBuf::from_str("tests/data/test_pops/populations.tsv")?;
    let chrom = "sq0";
    let begin = 0;
    let end = 1000;

    cmd.arg("loci")
        .arg(&input_file) // infile as a positional argument
        .arg(&output_prefix) // unique outfile as a positional argument
        .arg("-m")
        .arg("3")
        .arg("-M")
        .arg("15")
        .arg("-p")
        .arg(&pop_file);

    let output = cmd.output()?;
    debug!("Ran command");
    debug!("stdout: {}", String::from_utf8_lossy(&output.stdout));
    debug!("stderr: {}", String::from_utf8_lossy(&output.stderr));

    let success = d4_files_are_equal(expected_output_file.clone(), output_file.clone(), &chrom, begin, end)?;
    assert_eq!(true, success);
    fs::remove_file(&output_file)?;

    Ok(())
}

#[test]
fn loci_bgzf_test_no_pops_output_d4() -> Result<(), Box<dyn std::error::Error>> {
    //bamsim generate --num-bams 5 --num-chrs 1 --chr-length 1000 --min-depth 3 --max-depth 15 --min-mean-depth 6 --proportion 0.5 test_no_pops/
    let _ = env_logger::builder()
        .target(env_logger::Target::Stdout)
        .filter_level(log::LevelFilter::Trace)
        .is_test(true)
        .try_init();
    let mut cmd = Command::cargo_bin("clam")?;

    let input_file = PathBuf::from_str("tests/data/test_no_pops/merged.d4.gz")?;
    let output_prefix = PathBuf::from_str("tests/data/test_no_pops/output_d4")?;
    let output_file = output_prefix.with_extension("d4");
    let expected_output_file = PathBuf::from_str("tests/data/test_no_pops/truth_counts.d4")?;
    let chrom = "sq0";
    let begin = 0;
    let end = 1000;

    cmd.arg("loci")
        .arg(&input_file) // infile as a positional argument
        .arg(&output_prefix) // unique outfile as a positional argument
        .arg("-m")
        .arg("3")
        .arg("-M")
        .arg("15")
        .arg("-d")
        .arg("0.5")
        .arg("-u")
        .arg("6");

    let output = cmd.output()?;
    debug!("Ran command");
    debug!("stdout: {}", String::from_utf8_lossy(&output.stdout));
    debug!("stderr: {}", String::from_utf8_lossy(&output.stderr));

    let success = d4_files_are_equal(expected_output_file.clone(), output_file.clone(), &chrom, begin, end)?;
    assert_eq!(true, success);
    fs::remove_file(&output_file)?;

    Ok(())
}

#[test]
fn loci_bgzf_test_no_pops_output_bed() -> Result<(), Box<dyn std::error::Error>> {
    //bamsim generate --num-bams 5 --num-chrs 1 --chr-length 1000 --min-depth 3 --max-depth 15 --min-mean-depth 6 --proportion 0.5 test_no_pops/
    let _ = env_logger::builder()
        .target(env_logger::Target::Stdout)
        .filter_level(log::LevelFilter::Trace)
        .is_test(true)
        .try_init();
    let mut cmd = Command::cargo_bin("clam")?;

    let input_file = PathBuf::from_str("tests/data/test_no_pops/merged.d4.gz")?;
    let output_prefix = PathBuf::from_str("tests/data/test_no_pops/output_bed")?;
    let output_file = output_prefix.with_extension("bed");
    let expected_output_file = PathBuf::from_str("tests/data/test_no_pops/truth.bed")?;

    cmd.arg("loci")
        .arg(&input_file) // infile as a positional argument
        .arg(&output_prefix) // unique outfile as a positional argument
        .arg("-m")
        .arg("3")
        .arg("-M")
        .arg("15")
        .arg("-d")
        .arg("0.5")
        .arg("-u")
        .arg("6")
        .arg("--no-counts");

    let output = cmd.output()?;
    debug!("Ran command");
    debug!("stdout: {}", String::from_utf8_lossy(&output.stdout));
    debug!("stderr: {}", String::from_utf8_lossy(&output.stderr));

    let success = files_are_equal(expected_output_file.clone(), output_file.clone())?;
    assert_eq!(true, success);
    fs::remove_file(&output_file)?;

    Ok(())
}

#[test]
fn loci_test_no_pops_output_d4() -> Result<(), Box<dyn std::error::Error>> {
    //bamsim generate --num-bams 5 --num-chrs 1 --chr-length 1000 --min-depth 3 --max-depth 15 --min-mean-depth 6 --proportion 0.5 test_no_pops/
    let _ = env_logger::builder()
        .target(env_logger::Target::Stdout)
        .filter_level(log::LevelFilter::Trace)
        .is_test(true)
        .try_init();
    let mut cmd = Command::cargo_bin("clam")?;

    let input_file = PathBuf::from_str("tests/data/test_no_pops/merged.d4")?;
    let output_prefix = PathBuf::from_str("tests/data/test_no_pops/output_d4")?;
    let output_file = output_prefix.with_extension("d4");
    let expected_output_file = PathBuf::from_str("tests/data/test_no_pops/truth_counts.d4")?;
    let chrom = "sq0";
    let begin = 0;
    let end = 1000;

    cmd.arg("loci")
        .arg(&input_file) // infile as a positional argument
        .arg(&output_prefix) // unique outfile as a positional argument
        .arg("-m")
        .arg("3")
        .arg("-M")
        .arg("15")
        .arg("-d")
        .arg("0.5")
        .arg("-u")
        .arg("6");

    let output = cmd.output()?;
    debug!("Ran command");
    debug!("stdout: {}", String::from_utf8_lossy(&output.stdout));
    debug!("stderr: {}", String::from_utf8_lossy(&output.stderr));

    let success = d4_files_are_equal(expected_output_file.clone(), output_file.clone(), &chrom, begin, end)?;
    assert_eq!(true, success);
    fs::remove_file(&output_file)?;

    Ok(())
}

#[test]
fn loci_bgzf_test_pops() -> Result<(), Box<dyn std::error::Error>> {
    //./target/release/bamsim generate --num-bams 5 --num-chrs 1 --chr-length 1000 --min-depth 3 --max-depth 15 --min-mean-depth 6 --proportion 0.5 test_pops/pop1/bams/
    let _ = env_logger::builder()
        .target(env_logger::Target::Stdout)
        .filter_level(log::LevelFilter::Trace)
        .is_test(true)
        .try_init();
    let mut cmd = Command::cargo_bin("clam")?;

    let input_file = PathBuf::from_str("tests/data/test_pops/merged.d4.gz")?;
    let output_prefix = PathBuf::from_str("tests/data/test_pops/output_d4")?;
    let output_file = output_prefix.with_extension("d4");
    let expected_output_file = PathBuf::from_str("tests/data/test_pops/counts.d4")?;
    let pop_file = PathBuf::from_str("tests/data/test_pops/populations.tsv")?;
    let chrom = "sq0";
    let begin = 0;
    let end = 1000;

    cmd.arg("loci")
        .arg(&input_file) // infile as a positional argument
        .arg(&output_prefix) // unique outfile as a positional argument
        .arg("-m")
        .arg("3")
        .arg("-M")
        .arg("15")
        .arg("-p")
        .arg(&pop_file);

    let output = cmd.output()?;
    debug!("Ran command");
    debug!("stdout: {}", String::from_utf8_lossy(&output.stdout));
    debug!("stderr: {}", String::from_utf8_lossy(&output.stderr));

    let success = d4_files_are_equal(expected_output_file.clone(), output_file.clone(), &chrom, begin, end)?;
    assert_eq!(true, success);
    fs::remove_file(&output_file)?;

    Ok(())
}