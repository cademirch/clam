pub mod config;
pub mod roh;
pub mod utils;
pub mod vcf;
pub mod windows;
pub mod writer;

use crate::stat::config::StatConfig;
use crate::stat::windows::WindowStats;
use crate::stat::writer::Writer;
use color_eyre::Result;
use rayon::prelude::*;
use std::path::Path;


pub fn run_stat(config: StatConfig, outdir: &Path) -> Result<()> {
    let has_roh = config.roh_bed_path.is_some();
    let dxy = config.pop_map.num_populations() > 1;

    let mut results: Vec<(usize, Vec<WindowStats>)> = config
        .chunks
        .par_iter()
        .map(|(chunk)| {
            let query = config.create_query(chunk)?;
            let stats = query.process()?;
            Ok::<_, color_eyre::Report>((chunk.chunk_idx as usize, stats))
        })
        .collect::<Result<Vec<_>>>()?;

    
    results.sort_unstable_by_key(|(i, _)| *i);

    let mut writers = Writer::new(outdir, has_roh, dxy)?;
    for (_, stats) in results {
        for ws in stats {
            writers.write(ws)?;
        }
    }
    writers.flush()?;

    Ok(())
}