use crate::core::contig::ContigSet;
use color_eyre::{eyre::eyre, Result};
use indexmap::IndexMap;
use ndarray::Array2;
use std::marker::PhantomData;
use std::num::NonZeroUsize;
use std::ops::Range;
use std::path::{Path, PathBuf};
use std::sync::Arc;
use zarrs::array::codec::{BloscCodec, BloscCompressionLevel, BloscCompressor, BloscShuffleMode};

use zarrs::array::ElementOwned;
use zarrs::array::{Array, ArrayBuilder, DataType, FillValue};
use zarrs::array_subset::ArraySubset;
use zarrs::filesystem::FilesystemStore;
use zarrs::group::{Group, GroupBuilder};
use zarrs::storage::ReadableWritableListableStorage;

/// Multi-chromosome zarr storage with typed values
pub struct ChromosomeArrays<T: ElementOwned> {
    path: PathBuf,
    store: ReadableWritableListableStorage,
    contigs: ContigSet,
    column_names: Vec<String>,
    chunk_size: u64,
    _marker: PhantomData<T>,
}

impl<T: ElementOwned> ChromosomeArrays<T> {
    pub fn create(
        path: impl AsRef<Path>,
        contigs: ContigSet,
        column_names: Vec<String>,
        chunk_size: u64,
        data_type: DataType,
        fill_value: FillValue,
    ) -> Result<Self> {
        let path = path.as_ref().to_path_buf();
        let store = Self::open_store(&path)?;

        // Create root group and metadata
        let mut group = GroupBuilder::new().build(store.clone(), "/")?;
        let metadata = serde_json::json!({
            "contigs": contigs.iter()
                .map(|(name, len)| serde_json::json!({
                    "name": name,
                    "length": len
                }))
                .collect::<Vec<_>>(),
            "column_names": column_names,
            "chunk_size": chunk_size,
        });
        group
            .attributes_mut()
            .insert("clam_metadata".to_string(), metadata);
        group.store_metadata()?;

        // Derive typesize from data type for blosc codec
        let typesize = match &data_type {
            DataType::UInt8 | DataType::Int8 => 1,
            DataType::UInt16 | DataType::Int16 => 2,
            DataType::UInt32 | DataType::Int32 | DataType::Float32 => 4,
            DataType::UInt64 | DataType::Int64 | DataType::Float64 => 8,
            _ => 4, // default
        };

        // Create array for each chromosome
        for (chrom_name, chrom_length) in contigs.iter() {
            let mut array = ArrayBuilder::new(
                vec![chrom_length as u64, column_names.len() as u64],
                vec![chunk_size, column_names.len() as u64],
                data_type.clone(),
                fill_value.clone(),
            )
            .bytes_to_bytes_codecs(vec![Arc::new(BloscCodec::new(
                BloscCompressor::Zstd,
                BloscCompressionLevel::try_from(5).unwrap(),
                None, // blocksize (auto)
                BloscShuffleMode::Shuffle,
                Some(typesize), // typesize
            )?)])
            .build(store.clone(), &format!("/{}", chrom_name))?;

            // Add contig-specific attributes
            let contig_metadata = serde_json::json!({
                "contig": chrom_name,
                "length": chrom_length,
            });
            array
                .attributes_mut()
                .insert("contig_info".to_string(), contig_metadata);

            array.store_metadata()?;
        }

        Ok(Self {
            path,
            store,
            contigs,
            column_names,
            chunk_size,
            _marker: PhantomData,
        })
    }

    /// Open existing zarr
    pub fn open(path: impl AsRef<Path>) -> Result<Self> {
        let path = path.as_ref().to_path_buf();
        let store = Self::open_store(&path)?;

        // Read root metadata
        let group = Group::open(store.clone(), "/")?;
        let metadata = group
            .attributes()
            .get("clam_metadata")
            .ok_or_else(|| eyre!("Missing clam_metadata in zarr root"))?;

        let contigs_data = metadata["contigs"]
            .as_array()
            .ok_or_else(|| eyre!("Invalid contigs in metadata"))?;
        let contigs: Vec<_> = contigs_data
            .iter()
            .map(|c| {
                let name = c["name"].as_str().unwrap().to_string();
                let length = c["length"].as_u64().unwrap() as u32;
                crate::core::contig::Contig::new(name, length as usize)
            })
            .collect();
        let contigs = ContigSet::new(contigs);

        let column_names: Vec<String> = metadata["column_names"]
            .as_array()
            .ok_or_else(|| eyre!("Invalid column_names in metadata"))?
            .iter()
            .map(|v| v.as_str().unwrap().to_string())
            .collect();

        let chunk_size = metadata["chunk_size"]
            .as_u64()
            .ok_or_else(|| eyre!("Invalid chunk_size in metadata"))?;

        Ok(Self {
            path,
            store,
            contigs,
            column_names,
            chunk_size,
            _marker: PhantomData,
        })
    }

    fn open_store(path: &Path) -> Result<ReadableWritableListableStorage> {
        let store = Arc::new(FilesystemStore::new(path)?);
        Ok(store)
    }

    pub fn read_chunk(&self, chrom: &str, chunk_idx: u64) -> Result<Array2<T>> {
        let array = Array::open(self.store.clone(), &format!("/{}", chrom))?;
        let data = array.retrieve_chunk_ndarray(&[chunk_idx, 0])?;
        Ok(data.into_dimensionality()?)
    }

    /// Write a complete chunk
    pub fn write_chunk(&self, chrom: &str, chunk_idx: u64, data: Array2<T>) -> Result<()> {
        let array = Array::open(self.store.clone(), &format!("/{}", chrom))?;
        array.store_chunk_ndarray(&[chunk_idx, 0], data)?;
        Ok(())
    }

    pub fn contigs(&self) -> &ContigSet {
        &self.contigs
    }

    pub fn column_names(&self) -> &[String] {
        &self.column_names
    }

    pub fn chunk_size(&self) -> u64 {
        self.chunk_size
    }

    pub fn path(&self) -> &Path {
        &self.path
    }

    // Chunk utilities

    /// Convert position to chunk index
    pub fn position_to_chunk(&self, position: u32) -> u64 {
        (position as u64) / self.chunk_size
    }

    /// Get the start and end positions of a chunk
    pub fn chunk_bounds(&self, chunk_idx: u64, chrom_length: u32) -> (u32, u32) {
        let start = (chunk_idx * self.chunk_size) as u32;
        let end = ((chunk_idx + 1) * self.chunk_size).min(chrom_length as u64) as u32;
        (start, end)
    }

    /// Get the range of chunks that overlap a region
    pub fn overlapping_chunks(&self, start: u32, end: u32) -> Range<u64> {
        let start_chunk = self.position_to_chunk(start);
        let end_chunk = self.position_to_chunk(end.saturating_sub(1)) + 1;
        start_chunk..end_chunk
    }
}

pub type DepthArrays = ChromosomeArrays<u32>;
pub type CallableArrays = ChromosomeArrays<u8>;

impl DepthArrays {
    pub fn create_new(
        path: impl AsRef<Path>,
        contigs: ContigSet,
        sample_names: Vec<String>,
        chunk_size: u64,
    ) -> Result<Self> {
        Self::create(
            path,
            contigs,
            sample_names,
            chunk_size,
            DataType::UInt32,
            FillValue::from(0u32),
        )
    }
}

impl CallableArrays {
    pub fn create_new(
        path: impl AsRef<Path>,
        contigs: ContigSet,
        population_names: Vec<String>,
        chunk_size: u64,
    ) -> Result<Self> {
        Self::create(
            path,
            contigs,
            population_names,
            chunk_size,
            DataType::UInt8,
            FillValue::from(0u8),
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core::contig::Contig;
    use ndarray::Array2;
    use tempfile::TempDir;

    fn test_contigs() -> ContigSet {
        ContigSet::new(vec![
            Contig::new("chr1".to_string(), 1000),
            Contig::new("chr2".to_string(), 500),
        ])
    }

    #[test]
    fn test_create_and_open() {
        let dir = TempDir::new().unwrap();
        let path = dir.path().join("test.zarr");

        let contigs = test_contigs();
        let columns = vec!["sample1".to_string(), "sample2".to_string()];

        // Create
        DepthArrays::create_new(&path, contigs.clone(), columns.clone(), 100).unwrap();
        assert!(path.exists());

        // Open
        let arrays = DepthArrays::open(&path).unwrap();
        assert_eq!(arrays.contigs().len(), 2);
        assert_eq!(arrays.column_names(), &columns);
        assert_eq!(arrays.chunk_size(), 100);
    }

    #[test]
    fn test_write_and_read_chunk() {
        let dir = TempDir::new().unwrap();
        let path = dir.path().join("test.zarr");

        let arrays = DepthArrays::create_new(
            &path,
            test_contigs(),
            vec!["sample1".to_string(), "sample2".to_string()],
            100,
        )
        .unwrap();

        // Create test data (100 rows x 2 columns)
        let mut data = Array2::<u32>::zeros((100, 2));
        for i in 0..100 {
            data[[i, 0]] = 10;
            data[[i, 1]] = 20;
        }

        // Write and read back
        arrays.write_chunk("chr1", 0, data).unwrap();
        let read_data = arrays.read_chunk("chr1", 0).unwrap();

        assert_eq!(read_data.shape(), &[100, 2]);
        assert_eq!(read_data[[0, 0]], 10);
        assert_eq!(read_data[[99, 1]], 20);
    }

    #[test]
    fn test_chunk_utilities() {
        let dir = TempDir::new().unwrap();
        let path = dir.path().join("test.zarr");

        let arrays =
            DepthArrays::create_new(&path, test_contigs(), vec!["sample1".to_string()], 100)
                .unwrap();

        assert_eq!(arrays.position_to_chunk(0), 0);
        assert_eq!(arrays.position_to_chunk(99), 0);
        assert_eq!(arrays.position_to_chunk(100), 1);

        assert_eq!(arrays.chunk_bounds(0, 1000), (0, 100));
        assert_eq!(arrays.chunk_bounds(9, 1000), (900, 1000));

        assert_eq!(arrays.overlapping_chunks(0, 100), 0..1);
        assert_eq!(arrays.overlapping_chunks(50, 250), 0..3);
    }
}
