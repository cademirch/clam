use color_eyre::eyre::{bail, Context, Ok, Result};
use noodles::bgzf::Reader;
use noodles::vcf::io::indexed_reader;
use noodles::vcf::io::IndexedReader;
use noodles::vcf::Header;
use std::ffi::{OsStr, OsString};
use std::fs::File;
use std::path::{Path, PathBuf};
pub mod query;
pub mod variants;

/// Returns a path with a new dotted extension component appended to the end.
/// Note: does not check if the path is a file or directory; you should do that.
/// Source: https://internals.rust-lang.org/t/pathbuf-has-set-extension-but-no-add-extension-cannot-cleanly-turn-tar-to-tar-gz/14187/10

pub fn append_ext(ext: impl AsRef<OsStr>, path: PathBuf) -> PathBuf {
    let mut os_string: OsString = path.into();
    os_string.push(".");
    os_string.push(ext.as_ref());
    os_string.into()
}

pub fn build_vcf_reader(path: impl AsRef<Path>) -> Result<(IndexedReader<Reader<File>>, Header)> {
    let path = path.as_ref();
    let tbi_exists = append_ext(OsStr::new("tbi"), path.to_path_buf()).exists();
    let csi_exists = append_ext(OsStr::new("csi"), path.to_path_buf()).exists();

    if !tbi_exists && !csi_exists {
        bail!(
            "No index file found for {}. Expected .tbi or .csi",
            path.display()
        );
    }

    let mut reader = indexed_reader::Builder::default()
        .build_from_path(path)
        .wrap_err_with(|| format!("Failed to read VCF file: {}", path.display()))?;

    let header = reader.read_header().wrap_err("Failed to read VCF header")?;

    Ok((reader, header))
}
