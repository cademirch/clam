use super::DepthSource;
use crate::core::contig::{Contig, ContigSet};
use color_eyre::{eyre::WrapErr, Result, eyre::eyre};
use d4::index::D4IndexCollection;
use d4::ssio::D4TrackReader;
use log::warn;
use noodles::bgzf::{self, IndexedReader};
use std::fs::File;
use std::path::{Path, PathBuf};
use std::sync::OnceLock;
use std::sync::Mutex;
use std::collections::HashSet;

static WARNED_FILES: OnceLock<Mutex<HashSet<PathBuf>>> = OnceLock::new();

fn warn_missing_index(path: &Path) {
    let warned = WARNED_FILES.get_or_init(|| Mutex::new(HashSet::new()));
    let mut set = warned.lock().unwrap();
    if set.insert(path.to_path_buf()) {
        warn!(
            "SFI index not found for D4 file: {}, this will result in slower performance. \
            You can create the index by running: d4tools index build {}",
            path.display(),
            path.display()
        );
    }
}


pub struct D4Reader {
    inner: D4TrackReader<File>,
    sample_name: String,
    
}

impl D4Reader {
    /// Create a new D4Reader for a specific track
    ///
    /// # Arguments
    /// * `src` - Path to the D4 file
    /// * `sample_name` - The biological sample name to associate with this reader.
    ///                   This is what will be returned by `sample_name()` and used
    ///                   for downstream analysis.
    /// * `track_name` - The internal track path within the D4 file (e.g., "/sample1").
    ///                  If None, reads the default/only track in the file.
    ///                  In multitrack D4 files, we assume each track one sample.
    ///
    /// # Note
    /// For single-sample D4 files, `track_name` should be `None`.
    /// For multisample D4 files, `track_name` identifies which track to read,
    /// and `sample_name` is typically derived from the track name.
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
        let ic = D4IndexCollection::from_root_container(track_root);
        if ic.is_err() {
            warn_missing_index(&src_path);
        }

        Ok(Self {
            inner: d4_reader,
            sample_name: sample_name.to_string(),
            
        })
    }
    pub fn list_tracks<P: AsRef<Path>>(src: P) -> Result<Vec<PathBuf>> {
        let src_path = src.as_ref();
        let mut file = File::open(src_path)
            .wrap_err_with(|| format!("Failed to open D4 file: {}", src_path.display()))?;
        
        let mut tracks = Vec::new();
        d4::find_tracks(&mut file, |_| true, &mut tracks)
            .wrap_err("Failed to enumerate tracks in D4 file")?;
        
        Ok(tracks)
    }
    
    
    /// Create one reader per track in a multisample D4 file
    ///
    /// Each track in the D4 file becomes a separate DepthSource with the track name
    /// used as both the internal track identifier and the sample name.
    ///
    /// # Example
    /// A multisample D4 file with tracks "/sample1" and "/sample2" will return
    /// two readers where `sample_name()` returns "sample1" and "sample2" respectively.
    pub fn from_multisample_file<P: AsRef<Path>>(
        src: P,
    ) -> Result<Vec<Box<dyn DepthSource>>> {
        let src_path = src.as_ref();
        let track_paths = Self::list_tracks(src_path)?;
        
        
        track_paths
            .into_iter()
            .map(|track_path| {
                // Convert PathBuf to &str for track name
                let track_name = track_path
                    .to_str()
                    .ok_or_else(|| eyre!("Invalid track path: {:?}", track_path))?;
                
                // Use track name as sample name (e.g., "/sample1" -> "/sample1")
                // You could strip the leading "/" if desired: track_name.trim_start_matches('/')
                let sample_name = track_name;
                
                let reader = Self::new(src_path, sample_name, Some(track_name))?;
                Ok(Box::new(reader) as Box<dyn DepthSource>)
            })
            .collect()
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
            .map(|chrom| Contig::new(chrom.name.clone(), chrom.size))
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
        let ic = D4IndexCollection::from_root_container(track_root);
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
            .map(|chrom| Contig::new(chrom.name.clone(), chrom.size))
            .collect();

        ContigSet::new(contigs)
    }

    fn sample_name(&self) -> &str {
        &self.sample_name
    }
}
