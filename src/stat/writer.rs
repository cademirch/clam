use color_eyre::eyre::Ok;
use color_eyre::Result;
use csv;
use std::fs::File;
use std::path::Path;

use crate::stat::windows::WindowStats;

fn create_tsv_writer(outdir: &Path, filename: &str) -> Result<csv::Writer<File>> {
    let path = outdir.join(filename);
    let file = File::create(&path)?;
    Ok(csv::WriterBuilder::new().delimiter(b'\t').from_writer(file))
}

pub struct Writer {
    pub pi: csv::Writer<File>,
    pub dxy: Option<csv::Writer<File>>,
    pub fst: Option<csv::Writer<File>>,
    pub heterozygosity: Option<csv::Writer<File>>,
}

impl Writer {
    pub fn new<P: AsRef<Path>>(outdir: P, has_callable: bool, dxy: bool) -> Result<Self> {
        let pi = create_tsv_writer(outdir.as_ref(), "pi.tsv")?;
        let heterozygosity = if has_callable {
            Some(create_tsv_writer(outdir.as_ref(), "heterozygosity.tsv")?)
        } else {
            None
        };
        let (dxy, fst) = if dxy {
            (
                Some(create_tsv_writer(outdir.as_ref(), "dxy.tsv")?),
                Some(create_tsv_writer(outdir.as_ref(), "fst.tsv")?),
            )
        } else {
            (None, None)
        };
        Ok(Self {
            pi,
            dxy,
            fst,
            heterozygosity,
        })
    }

    pub fn write(&mut self, stats: WindowStats) -> Result<()> {
        for pi_stat in stats.pi {
            self.pi.serialize(pi_stat)?;
        }

        if let (Some(writer), Some(het_stats)) = (&mut self.heterozygosity, stats.heterozygosity) {
            for het_stat in het_stats {
                writer.serialize(het_stat)?;
            }
        }

        if let (Some(writer), Some(dxy)) = (&mut self.dxy, stats.dxy) {
            for dxy_stat in dxy {
                writer.serialize(dxy_stat)?;
            }
        }

        if let (Some(writer), Some(fst)) = (&mut self.fst, stats.fst) {
            for fst_stat in fst {
                writer.serialize(fst_stat)?;
            }
        }
        Ok(())
    }

    pub fn flush(&mut self) -> Result<()> {
        self.pi.flush()?;

        if let Some(writer) = &mut self.heterozygosity {
            writer.flush()?;
        }
        if let Some(writer) = &mut self.dxy {
            writer.flush()?;
        }
        if let Some(writer) = &mut self.fst {
            writer.flush()?;
        }
        Ok(())
    }
}
