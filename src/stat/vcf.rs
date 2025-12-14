use color_eyre::eyre::{Context, Ok, Result};
use noodles::bgzf::Reader;
use noodles::vcf::io::indexed_reader;
use noodles::vcf::io::IndexedReader;
use noodles::vcf::Header;
use std::fs::File;
use std::path::Path;

pub mod query;
pub mod variants;

pub fn build_vcf_reader(path: impl AsRef<Path>) -> Result<(IndexedReader<Reader<File>>, Header)> {
    let mut reader = indexed_reader::Builder::default()
        .build_from_path(path.as_ref())
        .wrap_err_with(|| format!("Failed to read VCF file: {}", path.as_ref().display()))?;

    let header = reader.read_header().wrap_err("Failed to read VCF header")?;

    Ok((reader, header))
}
