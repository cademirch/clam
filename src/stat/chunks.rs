use super::build_vcf_reader;
use anyhow::{Context, Result};
use bstr::BString;
use camino::Utf8PathBuf;
use fnv::FnvHashMap;
use noodles::core::{Position, Region};

const MIN_CHUNK_SIZE: usize = 100_000;

pub fn create_chunks(
    seqlens: FnvHashMap<String, usize>,
    vcf_path: Utf8PathBuf,
    window_size: usize,
) -> Result<Vec<(Region, bool)>> {
    let mut chunks = Vec::new();
    let mut current_chunk_size = 10_000_000;
    let (mut vcf_reader, header) = build_vcf_reader(vcf_path)?;

    for (chrom, length) in seqlens.iter() {
        let mut current_pos = 1;

        // First check for variants in the entire chromosome
        let begin = Position::new(current_pos).context("Invalid start position")?;
        let end = Position::new(*length).context("Invalid end position")?;
        let region = Region::new(BString::from(chrom.as_str()), begin..=end);
        let mut query = vcf_reader.query(&header, &region)?;

        // Get the position of the first variant, if any
        let first_variant_pos = if let Some(record) = query.next() {
            Some(record?.variant_start().unwrap().unwrap().get())
        } else {
            None
        };

        // If there is a first variant and it's not at the start, create windows up to it
        if let Some(variant_pos) = first_variant_pos {
            if variant_pos > 1 {
                let num_windows = (variant_pos - current_pos + window_size - 1) / window_size;
                for i in 0..num_windows {
                    let window_start = current_pos + (i * window_size);
                    if window_start >= variant_pos {
                        break;
                    }

                    let window_end = (window_start + window_size - 1).min(variant_pos - 1);
                    let begin = Position::new(window_start).context("Invalid start position")?;
                    let end = Position::new(window_end).context("Invalid end position")?;
                    let window_region = Region::new(BString::from(chrom.as_str()), begin..=end);

                    chunks.push((window_region, false));
                }
            }
            current_pos = variant_pos;
        }

        while current_pos <= *length {
            let begin = Position::new(current_pos).context("Invalid start position")?;
            let chunk_end = (current_pos + current_chunk_size - 1).min(*length);
            let end = Position::new(chunk_end).context("Invalid end position")?;

            let region = Region::new(BString::from(chrom.as_str()), begin..=end);
            let mut query = vcf_reader.query(&header, &region)?;

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

                    chunks.push((window_region, true));
                }

                current_chunk_size = (current_chunk_size / 2).max(MIN_CHUNK_SIZE);
            } else {
                let num_windows = (chunk_end - current_pos + window_size) / window_size;

                for i in 0..num_windows {
                    let window_start = current_pos + (i * window_size);
                    if window_start > chunk_end {
                        break;
                    }

                    let window_end = (window_start + window_size - 1).min(chunk_end);
                    let begin = Position::new(window_start).context("Invalid start position")?;
                    let end = Position::new(window_end).context("Invalid end position")?;
                    let window_region = Region::new(BString::from(chrom.as_str()), begin..=end);

                    chunks.push((window_region, false));
                }

                current_chunk_size = (current_chunk_size * 2).min(*length);
            }

            current_pos += current_chunk_size;
        }
    }

    Ok(chunks)
}
