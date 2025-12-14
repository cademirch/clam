pub mod vcf;
pub mod roh;
pub mod windows;
pub mod writer;
pub mod utils;
pub mod config;

use std::sync::mpsc::channel;
use rayon::prelude::*;
use std::path::Path;
use color_eyre::Result;
use crate::stat::writer::Writer;
use crate::stat::windows::WindowStats;
use crate::stat::config::StatConfig;


pub fn run_stat(config: StatConfig, outdir: &Path) -> Result<()> {
    
    let (tx, rx) = channel::<WindowStats>();
    
    let has_roh = config.roh_bed_path.is_some();
    let dxy = config.pop_map.num_populations() > 1;
    // Spawn dedicated writer thread
    let writer_handle = std::thread::spawn({
        let outdir = outdir.to_path_buf();
        move || -> Result<()> {
            let mut writers = Writer::new(&outdir, has_roh, dxy)?;
            
            // Receive and write stats as they arrive
            for window_stats in rx {
                writers.write(window_stats)?;
            }
            
            
            writers.flush()?;
            
            Ok(())
        }
    });
    
    
    let result = config.chunks
        .par_iter()
        .try_for_each(|chunk| {
            
            let query = config.create_query(chunk)?;
            
            
            let window_stats = query.process()?;
            
            
            for stats in window_stats {
                tx.send(stats)
                    .map_err(|_| color_eyre::eyre::eyre!("Failed to send window stats"))?;
            }
            
            Ok::<_, color_eyre::Report>(())
        });
    
    
    drop(tx);
    
    // Wait for writer thread to finish and handle any errors
    let writer_result = writer_handle.join()
        .map_err(|_| color_eyre::eyre::eyre!("Writer thread panicked"))?;
    
    
    result?;
    writer_result?;
    
    Ok(())
}