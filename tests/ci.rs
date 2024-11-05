use anyhow::Context;
use assert_cmd::prelude::*;
use log::debug;
use std::fs;
use std::io::Read;
use std::path::{Path, PathBuf};
use std::process::Command;
use std::str::FromStr;

#[test]
fn test_gzipped_d4() -> Result<(), Box<dyn std::error::Error>> {
    let _ = env_logger::builder()
        .target(env_logger::Target::Stdout)
        .filter_level(log::LevelFilter::Trace)
        .is_test(true)
        .try_init();
    let mut cmd = Command::cargo_bin("clam")?;

    // Generate a unique file name for the output
    let output_file = PathBuf::from_str("tests/data/test_gzipped.bed")?;

    // Specify the `loci` subcommand and add positional arguments for infile and outfile
    cmd.arg("loci")
        .arg("tests/data/merged.d4.gz") // infile as a positional argument
        .arg(&output_file) // unique outfile as a positional argument
        .arg("-m")
        .arg("3")
        .arg("-M")
        .arg("20")
        .arg("-d")
        .arg("0.66")
        .arg("-u")
        .arg("5");

    // Execute the command and capture the output
    let output = cmd.output()?;
    debug!("Ran command");
    debug!("stdout: {}", String::from_utf8_lossy(&output.stdout));
    debug!("stderr: {}", String::from_utf8_lossy(&output.stderr));
    // Read the generated output for comparison
    let mut generated_output = String::new();
    fs::File::open(&output_file)
        .context(format!("Couldn't open file: {}", output_file.display()))?
        .read_to_string(&mut generated_output)?;

    // Read the expected output from the truth file
    let expected_output_path = Path::new("tests/data/truth.bed");
    let mut expected_output = String::new();
    fs::File::open(expected_output_path)?.read_to_string(&mut expected_output)?;

    // Compare the generated output to the expected output
    assert_eq!(generated_output, expected_output);

    // Clean up by removing the generated output file
    fs::remove_file(&output_file)?;

    Ok(())
}

#[test]
fn test_gzipped_d4_pops() -> Result<(), Box<dyn std::error::Error>> {
    let _ = env_logger::builder()
        .target(env_logger::Target::Stdout)
        .filter_level(log::LevelFilter::Trace)
        .is_test(true)
        .try_init();
    let mut cmd = Command::cargo_bin("clam")?;

    // Generate a unique file name for the output
    let output_file = PathBuf::from_str("tests/data/test_gzipped_pops.d4")?;
    let pops_file = PathBuf::from_str("tests/data/populations.tsv")?;

    // Specify the `loci` subcommand and add positional arguments for infile and outfile
    cmd.arg("loci")
        .arg("tests/data/merged.d4.gz") // infile as a positional argument
        .arg(&output_file) // unique outfile as a positional argument
        .arg("-m")
        .arg("3")
        .arg("-M")
        .arg("20")
        .arg("-p")
        .arg(&pops_file);

    // Execute the command and capture the output
    let output = cmd.output()?;

    
    assert!(
        output.status.success(),
        "Command failed with stderr: {}",
        String::from_utf8_lossy(&output.stderr)
    );

    debug!("Ran command");
    debug!("stdout: {}", String::from_utf8_lossy(&output.stdout));
    debug!("stderr: {}", String::from_utf8_lossy(&output.stderr));

    
    let metadata = fs::metadata(&output_file)
        .context(format!("Couldn't get metadata for file: {}", output_file.display()))?;
    assert!(metadata.len() > 0, "Output file should not be empty.");

    
    fs::remove_file(&output_file)?;

    Ok(())
}

#[test]
fn test_d4() -> Result<(), Box<dyn std::error::Error>> {
    let _ = env_logger::builder()
        .target(env_logger::Target::Stdout)
        .filter_level(log::LevelFilter::Trace)
        .is_test(true)
        .try_init();
    let mut cmd = Command::cargo_bin("clam")?;

    // Generate a unique file name for the output
    let output_file = PathBuf::from_str("tests/data/test.bed")?;

    // Specify the `loci` subcommand and add positional arguments for infile and outfile
    cmd.arg("loci")
        .arg("tests/data/merged.d4") // infile as a positional argument
        .arg(&output_file) // unique outfile as a positional argument
        .arg("-m")
        .arg("3")
        .arg("-M")
        .arg("20")
        .arg("-d")
        .arg("0.66")
        .arg("-u")
        .arg("5");

    // Execute the command and capture the output
    let output = cmd.output()?;
    debug!("Ran command");
    debug!("stdout: {}", String::from_utf8_lossy(&output.stdout));
    debug!("stderr: {}", String::from_utf8_lossy(&output.stderr));
    // Read the generated output for comparison
    let mut generated_output = String::new();
    fs::File::open(&output_file)
        .context(format!("Couldn't open file: {}", output_file.display()))?
        .read_to_string(&mut generated_output)?;

    // Read the expected output from the truth file
    let expected_output_path = Path::new("tests/data/truth.bed");
    let mut expected_output = String::new();
    fs::File::open(expected_output_path)?.read_to_string(&mut expected_output)?;

    // Compare the generated output to the expected output
    assert_eq!(generated_output, expected_output);

    // Clean up by removing the generated output file
    fs::remove_file(&output_file)?;

    Ok(())
}