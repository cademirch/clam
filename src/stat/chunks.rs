use super::build_vcf_reader;
use anyhow::{Context, Result};
use bstr::BString;
use camino::Utf8PathBuf;
use noodles::core::{Position, Region};
use fnv::FnvHashMap;

const MIN_CHUNK_SIZE: usize = 100_000;

pub fn create_chunks(
    seqlens: FnvHashMap<String, usize>,
    vcf_path: Utf8PathBuf,
    window_size: usize,
) -> Result<Vec<Region>> {
    let mut chunks = Vec::new();
    let mut current_pos = 1;
    let mut current_chunk_size = 10_000_000;
    let (mut vcf_reader, header) = build_vcf_reader(vcf_path)?;

    for (chrom, length) in seqlens.iter() {
        let begin = Position::new(current_pos).context("Invalid start position")?;
        let end = Position::new(*length).context("Invalid end position")?;
        let region = Region::new(BString::from(chrom.as_str()), begin..=end);
        let mut query = vcf_reader.query(&header, &region)?;
        if let Some(record) = query.next() {
            current_pos = record?.variant_start().unwrap().unwrap().get();
        }
        while current_pos <= *length {
            let begin = Position::new(current_pos).context("Invalid start position")?;
            let chunk_end = (current_pos + current_chunk_size - 1).min(*length);
            let end = Position::new(chunk_end).context("Invalid end position")?;

            let region = Region::new(BString::from(chrom.as_str()), begin..=end);

            let mut query = vcf_reader.query(&header, &region)?;
            // let mut windows = Vec::new();

            if let Some(record) = query.next() {
                let variant_pos = record?.variant_start().unwrap().unwrap().get();
                let num_windows = (chunk_end - variant_pos + window_size) / window_size;

                for i in 0..num_windows {
                    let window_start = variant_pos + (i * window_size);
                    if window_start > chunk_end {
                        break;
                    }

                    let window_end = (window_start + window_size - 1).min(chunk_end);
                    let begin = Position::new(window_start).context("Invalid start position")?;
                    let end = Position::new(window_end).context("Invalid end position")?;
                    let window_region = Region::new(BString::from(chrom.as_str()), begin..=end);

                    chunks.push(window_region);
                }

                // Found variants - reduce chunk size to focus on dense regions
                current_chunk_size = (current_chunk_size / 2).max(MIN_CHUNK_SIZE);
            } else {
                // No variants - increase chunk size to skip sparse regions faster
                current_chunk_size = (current_chunk_size * 2).min(*length);
            }

            current_pos += current_chunk_size;
        }
    }

    Ok(chunks)
}
