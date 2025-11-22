use super::DepthSource;
use crate::core::contig::{Contig, ContigSet};
use color_eyre::{eyre::WrapErr, Result};
use d4::index::D4IndexCollection;
use d4::ssio::D4TrackReader;
use log::warn;
use noodles::bgzf::{self, IndexedReader};
use std::fs::File;
use std::path::{Path, PathBuf};

pub struct D4Reader {
    inner: D4TrackReader<File>,
    sample_name: String,
    src: PathBuf,
}

impl D4Reader {
    pub fn new<P: AsRef<Path>>(
        src: P,
        sample_name: &str,
        track_name: Option<&str>,
    ) -> Result<Self> {
        let src_path = src.as_ref().to_path_buf();

        let file = File::open(&src_path)
            .wrap_err_with(|| format!("Failed to open D4 file: {}", src_path.display()))?;

        let d4_reader = D4TrackReader::from_reader(file, track_name).map_err(|e| {
            color_eyre::eyre::eyre!(
                "Failed to create D4 track reader for file: {}, track: {:?}. Error: {}",
                src_path.display(),
                track_name,
                e
            )
        })?;

        // Check for SFI index
        let track_root = d4_reader.as_root();
        let ic = D4IndexCollection::from_root_container(&track_root);
        if ic.is_err() {
            warn!(
                "SFI index not found for D4 file: {}, this will result in slower performance. \
                You can create the index by running: d4tools index build {}",
                src_path.display(),
                src_path.display()
            );
        }

        Ok(Self {
            inner: d4_reader,
            sample_name: sample_name.to_string(),
            src: src_path,
        })
    }
}

impl DepthSource for D4Reader {
    fn read_depths(&mut self, chrom: &str, start: u32, end: u32) -> Result<Vec<u32>> {
        let mut view = self
            .inner
            .get_view(chrom, start, end)
            .wrap_err_with(|| format!("Failed to get view for {}:{}-{}", chrom, start, end))?;

        let len = (end - start) as usize;
        let mut depths = Vec::with_capacity(len);

        for pos in start..end {
            let (reported_pos, value) = view
                .next()
                .ok_or_else(|| {
                    color_eyre::eyre::eyre!("Unexpected end of D4 view at position {}", pos)
                })?
                .wrap_err("Failed to read D4 value")?;

            assert_eq!(
                reported_pos, pos,
                "D4 position mismatch: expected {}, got {}",
                pos, reported_pos
            );

            let depth: u32 = value
                .try_into()
                .wrap_err("Depth value doesn't fit in u32")?;

            depths.push(depth);
        }

        Ok(depths)
    }

    fn contigs(&self) -> ContigSet {
        let contigs: Vec<Contig> = self
            .inner
            .get_header()
            .chrom_list()
            .iter()
            .map(|chrom| Contig::new(chrom.name.clone(), chrom.size as usize))
            .collect();

        ContigSet::new(contigs)
    }

    fn sample_name(&self) -> &str {
        &self.sample_name
    }
}

pub struct BgzfD4Reader {
    inner: D4TrackReader<IndexedReader<File>>,
    sample_name: String,
    src: PathBuf,
}

impl BgzfD4Reader {
    pub fn new<P: AsRef<Path>>(
        src: P,
        sample_name: &str,
        track_name: Option<&str>,
    ) -> Result<Self> {
        let src_path = src.as_ref().to_path_buf();

        let indexed_reader = bgzf::indexed_reader::Builder::default()
            .build_from_path(&src_path)
            .wrap_err_with(|| format!("Failed to open bgzipped D4 file: {}", src_path.display()))?;

        let d4_reader = D4TrackReader::from_reader(indexed_reader, track_name).map_err(|e| {
            color_eyre::eyre::eyre!(
                "Failed to create D4 track reader for file: {}, track: {:?}. Error: {}",
                src_path.display(),
                track_name,
                e
            )
        })?;

        // Check for SFI index
        let track_root = d4_reader.as_root();
        let ic = D4IndexCollection::from_root_container(&track_root);
        if ic.is_err() {
            warn!(
                "SFI index not found for D4 file: {}, this will result in slower performance. \
                You can create the index by running: d4tools index build {}",
                src_path.display(),
                src_path.display()
            );
        }

        Ok(Self {
            inner: d4_reader,
            sample_name: sample_name.to_string(),
            src: src_path,
        })
    }
}

impl DepthSource for BgzfD4Reader {
    fn read_depths(&mut self, chrom: &str, start: u32, end: u32) -> Result<Vec<u32>> {
        let mut view = self
            .inner
            .get_view(chrom, start, end)
            .wrap_err_with(|| format!("Failed to get view for {}:{}-{}", chrom, start, end))?;

        let len = (end - start) as usize;
        let mut depths = Vec::with_capacity(len);

        for pos in start..end {
            let (reported_pos, value) = view
                .next()
                .ok_or_else(|| {
                    color_eyre::eyre::eyre!("Unexpected end of D4 view at position {}", pos)
                })?
                .wrap_err("Failed to read D4 value")?;

            assert_eq!(
                reported_pos, pos,
                "D4 position mismatch: expected {}, got {}",
                pos, reported_pos
            );

            let depth: u32 = value
                .try_into()
                .wrap_err("Depth value doesn't fit in u32")?;

            depths.push(depth);
        }

        Ok(depths)
    }

    fn contigs(&self) -> ContigSet {
        let contigs: Vec<Contig> = self
            .inner
            .get_header()
            .chrom_list()
            .iter()
            .map(|chrom| Contig::new(chrom.name.clone(), chrom.size as usize))
            .collect();

        ContigSet::new(contigs)
    }

    fn sample_name(&self) -> &str {
        &self.sample_name
    }
}
