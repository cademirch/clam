use std::collections::HashSet;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

use camino::Utf8PathBuf;
use color_eyre::{eyre::bail, Result};

use noodles::bed;
use noodles::core::Region;

pub fn count_combinations(n: u32, r: u32) -> u32 {
    if r > n {
        0
    } else {
        (1..=r).fold(1, |acc, val| acc * (n - val + 1) / val)
    }
}

pub fn get_exclude_chromosomes(
    exclude: &Option<Vec<String>>,
    exclude_file: &Option<Utf8PathBuf>,
) -> Result<Option<HashSet<String>>> {
    if let Some(file_path) = exclude_file {
        // Read chromosomes from the file and collect them into a HashSet
        let file = File::open(file_path)?;
        let reader = BufReader::new(file);
        let exclude_set: HashSet<String> = reader
            .lines()
            .filter_map(Result::ok) // Ignore any lines that fail to read
            .collect();
        Ok(Some(exclude_set))
    } else if let Some(chromosomes) = exclude {
        // Convert Vec<String> to HashSet<String> directly
        Ok(Some(chromosomes.iter().cloned().collect()))
    } else {
        // If neither option is provided, return an empty HashSet
        Ok(None)
    }
}

pub fn read_bed_regions<P: AsRef<Path>>(bed_file: P) -> Result<Vec<Region>> {
    let mut res = vec![];
    let mut rdr = bed::io::reader::Builder::<3>
        .build_from_path(bed_file.as_ref())
        .unwrap_or_else(|_| panic!("Failed to open bed file: {}",
            bed_file.as_ref().display()));
    let mut record = bed::Record::<3>::default();

    while match rdr.read_record(&mut record) {
        Ok(0) => false, // no more records break loop
        Ok(_) => true,  // continue loop
        Err(e) => {
            bail!(
                "Error reading record: {}. Offending record: {:?}",
                e,
                record
            );
        }
    } {
        let chrom = record.reference_sequence_name();
        let start = match record.feature_start() {
            Ok(pos) => pos,
            Err(e) => {
                bail!("Bad bed record start: {:?}. Error: {:?}", record, e);
            }
        };
        let end = match record.feature_end() {
            Some(Ok(pos)) => pos,
            Some(Err(e)) => bail!("Bad bed record end: {:?}. Error: {:?}", record, e),
            None => bail!("Bed record missing end: {:?}", record),
        };
        let region = Region::new(chrom, start..=end);

        res.push(region);
    }

    Ok(res)
}
