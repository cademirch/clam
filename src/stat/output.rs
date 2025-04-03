use anyhow::{Context, Result};
use csv;
use serde::Serialize;

use std::fs::File;
use std::path::Path;

pub fn create_output_file<P: AsRef<Path>>(outdir: P, file_name: &str) -> Result<csv::Writer<File>> {
    let output_path = outdir.as_ref().join(file_name);
    let output_file = File::create(&output_path).context(format!(
        "Couldn't create output file: {}",
        &output_path.display()
    ))?;

    let writer = csv::WriterBuilder::new()
        .delimiter(b'\t')
        .from_writer(output_file);

    Ok(writer)
}
#[derive(Debug, Serialize)]
pub struct HetRecord {
    pub chrom: String,
    pub start: u32,
    pub end: u32,
    pub callable_bases: u32,
    pub bases_in_roh:u32,
    pub sample_name: String,
    pub count_het_sites:u32,
}

impl HetRecord {
    pub fn new(
        chrom: &str,
        start: u32,
        end: u32,
        callable_bases: u32,
        bases_in_roh: u32,
        sample_name: &str,
        count_het_sites: u32,
    ) -> Self {
        HetRecord {
            chrom: chrom.to_string(),
            start,
            end,
            callable_bases,
            bases_in_roh,
            sample_name: sample_name.to_string(),
            count_het_sites,
        }
    }
}

#[derive(Serialize)]
pub struct PiRecord {
    population_name: String,
    chrom: String,
    start: u32,
    end: u32,
    pi: f32,
    comps: u32,
    diffs: u32,
}

impl PiRecord {
    pub fn new(
        population_name: &str,
        chrom: &str,
        start: u32,
        end: u32,
        pi: f32,
        comps: u32,
        diffs: u32,
    ) -> Self {
        let population_name = population_name.to_string();
        let chrom = chrom.to_string();
        Self {
            population_name,
            chrom,
            start,
            end,
            pi,
            comps,
            diffs,
        }
    }
}
#[derive(Serialize)]
pub struct DxyRecord {
    population1_name: String,
    population2_name: String,
    chrom: String,
    start: u32,
    end: u32,
    dxy: f32,
    comparisons: u32,
    differences: u32,
}

impl DxyRecord {
    pub fn new(
        population1_name: &str,
        population2_name: &str,
        chrom: &str,
        start: u32,
        end: u32,
        dxy: f32,
        comparisons: u32,
        differences: u32,
    ) -> Self {
        Self {
            population1_name: population1_name.to_string(),
            population2_name: population2_name.to_string(),
            chrom: chrom.to_string(),
            start,
            end,
            dxy,
            comparisons,
            differences,
        }
    }
}

#[derive(Serialize)]
pub struct FstRecord {
    population1_name: String,
    population2_name: String,
    chrom: String,
    start: u32,
    end: u32,
    fst: f32,
    
}

impl FstRecord {
    pub fn new(
        population1_name: &str,
        population2_name: &str,
        chrom: &str,
        start: u32,
        end: u32,
        fst: f32,
        
    ) -> Self {
        Self {
            population1_name: population1_name.to_string(),
            population2_name: population2_name.to_string(),
            chrom: chrom.to_string(),
            start,
            end,
            fst,
        }
    }
}