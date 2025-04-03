use std::fs::File;
use std::io::Write;
use std::path::{Path, PathBuf};

use super::regions::CallableRegion;
use anyhow::{Context, Result};
use d4::index::D4IndexCollection;
use d4::ptab::PTablePartitionWriter;
use d4::stab::SecondaryTablePartWriter;
use d4::{Chrom, D4FileBuilder, D4FileMerger, D4FileWriter, Dictionary};
use log::trace;
use rayon::prelude::*;

use tempfile::NamedTempFile;

pub fn write_d4_parallel<P: AsRef<Path>>(
    regions: &[(String, u32, Vec<CallableRegion>)],
    chroms: Vec<Chrom>,
    output_path: Option<P>,
) -> Result<PathBuf> {
    // Initialize the D4 file writer
    let output_path: PathBuf = if let Some(output_path) = output_path {
        output_path.as_ref().to_path_buf()
    } else {
        let temp_file = NamedTempFile::new()?;
        temp_file.into_temp_path().to_path_buf()
    };

    let mut builder = D4FileBuilder::new(output_path.clone());
    builder.append_chrom(chroms.into_iter());
    builder.set_denominator(1.0);
    builder.set_dictionary(Dictionary::new_simple_range_dict(0, 1)?);

    // Create the D4 file writer
    let mut d4_writer: D4FileWriter = builder.create()?;
    let partitions = d4_writer.parallel_parts(Some(10_000_000))?;
    trace!("Created {} partitions for writing", partitions.len());

    let partitions_result: Result<()> = partitions
        .into_par_iter()
        .map(|(mut partition, mut secondary_table)| {
            let (chrom, left, right) = partition.region();
            let chrom = chrom.to_string();
            let mut primary_encoder = partition.make_encoder();
            let mut last = left;

            // Encapsulate encoding logic in a closure
            let mut write_value = |pos: u32, value: i32| {
                if !primary_encoder.encode(pos as usize, value) {
                    secondary_table.encode(pos, value).unwrap();
                    // trace!("Secondary encode: Wrote {} at {}:{}", value, chrom, pos);
                } else {
                    // trace!("Primary encode: Wrote {} at {}:{}", value, chrom, pos);
                }
            };

            // Process regions within the current partition
            if let Some((_, begin_offset, region_list)) = regions
                .iter()
                .find(|(region_chrom, _, _)| *region_chrom == chrom)
            {
                for region in region_list {
                    let start = region.begin + *begin_offset;
                    let end = region.end + *begin_offset;

                    if end <= left || start >= right {
                        // Skip regions outside of the partition bounds
                        continue;
                    }

                    // Clip the region to fit within the partition bounds
                    let clipped_start = start.max(left);
                    let clipped_end = end.min(right);

                    for pos in last..clipped_start {
                        write_value(pos, 0);
                    }

                    for pos in clipped_start..clipped_end {
                        write_value(pos, region.count as i32);
                    }

                    last = clipped_end;
                }
            }

            // Fill remaining positions with zero
            for pos in last..right {
                write_value(pos, 0);
            }

            // Flush the secondary table
            secondary_table.flush()?;
            secondary_table.finish()?;
            Ok(())
        })
        .collect();

    // Check for any errors during the parallel operation
    if let Err(e) = partitions_result {
        return Err(e.into());
    }

    drop(d4_writer);
    log::info!("D4 file writing complete.");

    let mut index_collection = D4IndexCollection::open_for_write(output_path.clone())
        .with_context(|| {
            format!(
                "Failed to open D4 index collection for writing at {:?}",
                output_path
            )
        })?;

    index_collection
        .create_secondary_frame_index()
        .with_context(|| "Failed to create secondary frame index for the D4 index collection")?;
    Ok(output_path)
}

pub fn merge_d4_files<P: AsRef<Path>>(outpath: P, input: Vec<P>, names: Vec<&str>) -> Result<()> {
    let mut merger = D4FileMerger::new(outpath);
    for (file, name) in input.iter().zip(names.iter()) {
        merger = merger.add_input_with_tag(file, name)
    }
    merger.merge()?;
    Ok(())
}

pub fn write_bed<P: AsRef<Path>>(
    output_path: P,
    regions: &[(String, u32, Vec<CallableRegion>)],
) -> Result<()> {
    let mut file = File::create(output_path)?;

    for (chrom, begin, callable_regions) in regions {
        for region in callable_regions.into_iter() {
            writeln!(
                file,
                "{}\t{}\t{}",
                chrom,
                region.begin + begin,
                region.end + begin
            )?;
        }
    }

    Ok(())
}
