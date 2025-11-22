use bstr::ByteSlice;
use color_eyre::eyre::{Context, ContextCompat, OptionExt};
use color_eyre::{
    eyre::{bail, ensure, eyre, WrapErr},
    Result,
};
use flate2::read::{self, GzDecoder};
use noodles::bed;
use noodles::vcf::Header;
use rust_lapper::{Interval, Lapper};
use std::collections::HashMap;
use std::fs::File;
use std::io::{self, BufRead, BufReader, Read};
use std::path::Path;

type RohIndex = HashMap<String, Lapper<u32, usize>>;

fn open_bed<P: AsRef<Path>>(path: P) -> Result<Box<dyn BufRead>> {
    let file = File::open(&path).map_err(|e| match e.kind() {
        io::ErrorKind::NotFound => {
            eyre!("File not found: {}", path.as_ref().display())
        }
        io::ErrorKind::PermissionDenied => {
            eyre!("Permission denied: {}", path.as_ref().display())
        }
        _ => eyre!("Cannot open file: {} ({})", path.as_ref().display(), e),
    })?;
    let reader: Box<dyn BufRead> = if path
        .as_ref()
        .extension()
        .and_then(|ext| ext.to_str())
        .map_or(false, |ext| ext == "gz")
    {
        Box::new(BufReader::new(GzDecoder::new(file)))
    } else {
        Box::new(BufReader::new(file))
    };

    Ok(reader)
}

fn create_roh_index<R: BufRead>(bed_file: R, vcf_header: &Header) -> Result<RohIndex> {
    let mut bed_reader = bed::io::reader::Reader::<4, _>::new(bed_file);
    let contigs = vcf_header.contigs();
    let samples = vcf_header.sample_names();

    let mut contig_intervals: HashMap<String, Vec<Interval<u32, usize>>> = HashMap::new();
    let mut record = bed::Record::default();

    while bed_reader
        .read_record(&mut record)
        .wrap_err_with(|| format!("Failed to read bed record: {:?}", record))?
        != 0
    {
        let chrom = record.reference_sequence_name().to_str_lossy().into_owned();

        let contig_key = if contigs.contains_key(&chrom) {
            chrom
        } else {
            bail!("BED contig '{}' not found in VCF header", chrom);
        };

        let start = record
            .feature_start()
            .wrap_err_with(|| format!("Invalid start in bed record: {:?}", record))?;

        let end = record
            .feature_end()
            .transpose()
            .wrap_err_with(|| format!("Invalid end in bed record: {:?}", record))?
            .ok_or_else(|| eyre!("Missing end in bed record: {:?}", record))?;

        let sample_name = record
            .name()
            .ok_or_else(|| eyre!("Missing sample name in bed record: {:?}", record))?
            .to_str_lossy()
            .into_owned();

        let sample_idx = samples.get_index_of(&sample_name).ok_or_eyre(format!(
            "BED sample: {} not found in VCF header.",
            sample_name
        ))?;

        contig_intervals
            .entry(contig_key)
            .or_default()
            .push(Interval {
                // noodles already converted BED's 0-based coords to 1-based Position
                // e.g., BED [0, 100) becomes Position(1) to Position(100)
                start: start.get() as u32,
                // Lapper uses half-open intervals, so Position [1, 100] inclusive
                // needs to become [1, 101) half-open to include position 100
                stop: end.get() as u32 + 1,
                val: sample_idx,
            });
    }

    let mut lappers = HashMap::new();
    for (contig, intervals) in contig_intervals {
        lappers.insert(contig, Lapper::new(intervals));
    }

    Ok(lappers)
}
