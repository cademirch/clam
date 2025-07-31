use std::fs::File;
use std::io::{BufReader, Read};
use std::path::{Path, PathBuf};
use std::process::Command;
use std::str::FromStr;

use anyhow::{Context, Result};
use assert_cmd::prelude::*;
use d4::find_tracks_in_file;
use d4::ssio::D4TrackReader;
use log::debug;

/// Compare two files byte-by-byte to check if they are the same.
fn files_are_equal<P: AsRef<Path> + Clone>(file1: P, file2: P) -> Result<bool> {
    let f1 = File::open(file1.clone())
        .context(format!("{} does not exist", file1.as_ref().display()))?;
    let f2 = File::open(file2.clone())
        .context(format!("{} does not exist", file2.as_ref().display()))?;

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

fn d4_files_are_equal<P: AsRef<Path> + Clone>(
    file1: P,
    file2: P,
    chrom: &str,
    begin: u32,
    end: u32,
) -> Result<bool> {
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
            f1_data[idx].push(value as u32); // Assuming values are within u32 range
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
    for (idx, (track1, track2)) in f1_data.iter().zip(f2_data.iter()).enumerate() {
        if track1 != track2 {
            for (pos, (value1, value2)) in track1.iter().zip(track2.iter()).enumerate() {
                if value1 != value2 {
                    debug!(
                        "Mismatch at track index {}, position {}: track1 value = {}, track2 value = {}",
                        idx, pos, value1, value2
                    );
                    return Ok(false); // Mismatch found
                }
            }
        }
    }

    Ok(true) // All data matched
}

fn init_logger() {
    let _ = env_logger::builder()
        .target(env_logger::Target::Stdout)
        .filter_level(log::LevelFilter::Trace)
        .is_test(true)
        .try_init();
}

fn file_paths(
    test_case: &str,
    outdir: &str,
    truth_file: &str,
    bgzip: bool,
) -> (PathBuf, PathBuf, PathBuf) {
    let input_file_str = if bgzip {
        format!("tests/data/loci/{}/filelist.txt", test_case)
    } else {
        format!("tests/data/loci/{}/merged.d4", test_case)
    };
    let input_file = PathBuf::from_str(&input_file_str).unwrap();
    let output_dir = PathBuf::from_str(outdir).unwrap();
    let expected_file =
        PathBuf::from_str(&format!("tests/data/loci/{}/{}", test_case, truth_file)).unwrap();
    (input_file, output_dir, expected_file)
}

#[test]
fn test_merged_pops() -> Result<()> {
    init_logger();
    let input_file = "tests/data/loci/test_pops/merged.d4";
    let output_dir = "tests/data/loci/test_pops/output";
    let extra_args = vec![
        "-p",
        "tests/data/loci/test_pops/populations.tsv",
        "--merged",
    ];

    std::fs::remove_dir_all(output_dir);
    let mut cmd = Command::cargo_bin("clam").unwrap();
    let clam_command = cmd
        .arg("loci")
        .arg(input_file)
        .arg("-o")
        .arg(output_dir)
        .arg("-m")
        .arg("3")
        .arg("-M")
        .arg("15")
        .args(&extra_args);

    let output = clam_command.output()?;
    log::debug!("stdout: {}", String::from_utf8_lossy(&output.stdout));
    log::debug!("stderr: {}", String::from_utf8_lossy(&output.stderr));

    let output_d4 = "tests/data/loci/test_pops/output/callable_sites.d4";
    let truth_d4 = "tests/data/loci/test_pops/counts.d4";

    let chrom = "sq0";
    let begin = 0;
    let end = 1000;
    assert!(d4_files_are_equal(truth_d4, output_d4, chrom, begin, end)?);
    std::fs::remove_dir_all(output_dir)?;
    Ok(())
}

#[test]
fn test_merged_no_pops() -> Result<()> {
    init_logger();
    let input_file = "tests/data/loci/test_no_pops/merged.d4";
    let output_dir = "tests/data/loci/test_no_pops/output";
    let extra_args = vec!["--merged"];
    std::fs::remove_dir_all(output_dir);
    let mut cmd = Command::cargo_bin("clam").unwrap();
    let clam_command = cmd
        .arg("loci")
        .arg(input_file)
        .arg("-o")
        .arg(output_dir)
        .arg("-m")
        .arg("3")
        .arg("-M")
        .arg("15")
        .args(&extra_args);

    let output = clam_command.output()?;
    log::debug!("stdout: {}", String::from_utf8_lossy(&output.stdout));
    log::debug!("stderr: {}", String::from_utf8_lossy(&output.stderr));

    let output_d4 = "tests/data/loci/test_no_pops/output/callable_sites.d4";
    let truth_d4 = "tests/data/loci/test_no_pops/truth_counts.d4";

    let chrom = "sq0";
    let begin = 0;
    let end = 1000;
    assert!(d4_files_are_equal(truth_d4, output_d4, chrom, begin, end)?);
    std::fs::remove_dir_all(output_dir)?;
    Ok(())
}

#[test]
fn test_multid4_no_pops() -> Result<()> {
    init_logger();
    let input_file = "tests/data/loci/test_no_pops/bgzf/filelist.txt";
    let output_dir = "tests/data/loci/test_no_pops/bgzf/output";

    std::fs::remove_dir_all(output_dir);
    let mut cmd = Command::cargo_bin("clam").unwrap();
    let clam_command = cmd
        .arg("loci")
        .arg("-f")
        .arg(input_file)
        .arg("-o")
        .arg(output_dir)
        .arg("-m")
        .arg("3")
        .arg("-M")
        .arg("15");

    let output = clam_command.output()?;
    log::debug!("stdout: {}", String::from_utf8_lossy(&output.stdout));
    log::debug!("stderr: {}", String::from_utf8_lossy(&output.stderr));

    let output_d4 = "tests/data/loci/test_no_pops/bgzf/output/callable_sites.d4";
    let truth_d4 = "tests/data/loci/test_no_pops/bgzf/truth_counts.d4";
    let output_bed = "tests/data/loci/test_no_pops/bgzf/output/callable_sites.bed";
    let chrom = "sq0";
    let begin = 0;
    let end = 1000;
    assert!(d4_files_are_equal(truth_d4, output_d4, chrom, begin, end)?);
    assert_ne!(PathBuf::from_str(output_bed).unwrap().exists(), true);
    std::fs::remove_dir_all(output_dir)?;
    Ok(())
}
