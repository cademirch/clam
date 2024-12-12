use super::*;
use anyhow::{bail, Context, Result};
use fnv::FnvHashMap;
use noodles::vcf::{
    self,
    variant::record::samples::{keys::key, series::Value, Series},
};

/// Count alleles in a single vcf record
pub fn count_alleles(
    record: &vcf::Record,
    header: &vcf::Header,
    sample_to_pop: &FnvHashMap<usize, usize>,
    counts: &mut Vec<[u32; 2]>,
    samples_in_roh: HashSet<usize>,
) -> Result<()> {
    let samples = record.samples();
    let gt_series = samples
        .select(key::GENOTYPE)
        .ok_or_else(|| anyhow!("Malformed variant record: {:?}", record))?;

    for (sample_idx, value) in gt_series.iter(&header).enumerate() {
        let Some(Value::Genotype(genotype)) = value? else {
            bail!(
                "Invalid GT value for sample at index {}, in variant record: {:?}",
                sample_idx,
                record
            )
        };
        let pop_idx = sample_to_pop
            .get(&sample_idx)
            .ok_or_else(|| anyhow!("Unknown population for sample index {}", sample_idx))?;
        let mut ref_count = 0;
        let mut alt_count = 0;
        for result in genotype.iter() {
            let allele = result.with_context(|| {
                format!(
                    "Failed to parse alleles for sample at index {}, in variant record: {:?}",
                    sample_idx, record
                )
            })?;

            let (position, _) = allele;

            match position {
                Some(0) => ref_count += 1,
                Some(1) => alt_count += 1,
                None => (),
                _ => bail!(
                    "Unexpected allele value in sample index {}: {:?}",
                    sample_idx,
                    position
                ),
            }
        }

        if samples_in_roh.contains(&sample_idx) {
            trace!(
                "ROH sampidx: {}: ref:{} alt: {}",
                sample_idx,
                ref_count,
                alt_count
            );
            match ref_count.cmp(&alt_count) {
                std::cmp::Ordering::Greater => counts[*pop_idx][0] += 1,
                std::cmp::Ordering::Less => counts[*pop_idx][1] += 1,
                std::cmp::Ordering::Equal => (),
            }
        } else {
            counts[*pop_idx][0] += ref_count;
            counts[*pop_idx][1] += alt_count;
        }
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use anyhow::Result;
    use fnv::FnvHashMap;
    use noodles::vcf::io::Reader;
    use std::io::Cursor;

    #[test]
    fn test_count_alleles_single_population() -> Result<()> {
        let data: &str = concat!(
            "##fileformat=VCFv4.3\n",
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample1\tsample2\n",
            "sq0\t1\t.\tA\tG\t.\tPASS\t.\tGT\t0/0\t0/1\n"
        );

        let mut sample_to_pop = FnvHashMap::default();
        sample_to_pop.insert(0, 0); // sample1 -> pop0
        sample_to_pop.insert(1, 0); // sample2 -> pop0

        let data_cursor = Cursor::new(data);
        let mut reader = Reader::new(data_cursor);
        let actual_header = reader.read_header()?;

        for result in reader.records() {
            let record = result?;
            let mut counts: Vec<[u32; 2]> = vec![[0; 2]; 1]; // 1 population

            count_alleles(
                &record,
                &actual_header,
                &sample_to_pop,
                &mut counts,
                HashSet::new(),
            )?;

            assert_eq!(counts[0][0], 3); // 2 from sample1, 1 from sample2
            assert_eq!(counts[0][1], 1); // 1 from sample2
        }

        Ok(())
    }

    #[test]
    fn test_count_alleles_haploid_single_population() -> Result<()> {
        let data: &str = concat!(
            "##fileformat=VCFv4.3\n",
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample1\tsample2\n",
            "sq0\t1\t.\tA\tG\t.\tPASS\t.\tGT\t0\t1\n"
        );

        let mut sample_to_pop = FnvHashMap::default();
        sample_to_pop.insert(0, 0); // sample1 -> pop0
        sample_to_pop.insert(1, 0); // sample2 -> pop0

        let data_cursor = Cursor::new(data);
        let mut reader = Reader::new(data_cursor);
        let actual_header = reader.read_header()?;

        for result in reader.records() {
            let record = result?;
            let mut counts: Vec<[u32; 2]> = vec![[0; 2]; 1]; // 1 population

            count_alleles(
                &record,
                &actual_header,
                &sample_to_pop,
                &mut counts,
                HashSet::new(),
            )?;

            assert_eq!(counts[0][0], 1); // 1 from sample1
            assert_eq!(counts[0][1], 1); // 1 from sample2
        }

        Ok(())
    }

    #[test]
    fn test_count_alleles_two_populations() -> Result<()> {
        let data: &str = concat!(
            "##fileformat=VCFv4.3\n",
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample1\tsample2\n",
            "sq0\t1\t.\tA\tG\t.\tPASS\t.\tGT\t0/0\t0/1\n"
        );

        let mut sample_to_pop = FnvHashMap::default();
        sample_to_pop.insert(0, 0); // sample1 -> pop0
        sample_to_pop.insert(1, 1); // sample2 -> pop1

        let data_cursor = Cursor::new(data);
        let mut reader = Reader::new(data_cursor);
        let actual_header = reader.read_header()?;

        for result in reader.records() {
            let record = result?;
            let mut counts: Vec<[u32; 2]> = vec![[0; 2]; 2]; // 2 populations

            count_alleles(
                &record,
                &actual_header,
                &sample_to_pop,
                &mut counts,
                HashSet::new(),
            )?;

            // Population 0
            assert_eq!(counts[0][0], 2); // 2 from sample1
            assert_eq!(counts[0][1], 0); // 0 from sample1

            // Population 1
            assert_eq!(counts[1][0], 1); // 1 from sample2
            assert_eq!(counts[1][1], 1); // 1 from sample2
        }

        Ok(())
    }

    #[test]
    fn test_count_alleles_with_roh() -> Result<()> {
        let data: &str = concat!(
            "##fileformat=VCFv4.3\n",
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample1\tsample2\n",
            "sq0\t1\t.\tA\tG\t.\tPASS\t.\tGT\t0/0\t0/1\n",
            "sq0\t2\t.\tA\tG\t.\tPASS\t.\tGT\t1/1\t0/1\n"
        );

        let mut sample_to_pop = FnvHashMap::default();
        sample_to_pop.insert(0, 0); // sample1 -> pop0
        sample_to_pop.insert(1, 0); // sample2 -> pop0

        let roh_samples: HashSet<usize> = vec![0].into_iter().collect(); // Only sample1 is in ROH

        let data_cursor = Cursor::new(data);
        let mut reader = Reader::new(data_cursor);
        let actual_header = reader.read_header()?;

        for result in reader.records() {
            let record = result?;
            let mut counts: Vec<[u32; 2]> = vec![[0; 2]; 1]; // 1 population

            // Pass the counts array to count_alleles with ROH samples
            count_alleles(
                &record,
                &actual_header,
                &sample_to_pop,
                &mut counts,
                roh_samples.clone(),
            )?;

            // Assertions depend on the record and ROH status
            match record.variant_start().unwrap().unwrap().get() {
                1 => {
                    assert_eq!(counts[0][0], 2);
                    assert_eq!(counts[0][1], 1);
                }
                2 => {
                    assert_eq!(counts[0][0], 1);
                    assert_eq!(counts[0][1], 2);
                }
                _ => panic!("Unexpected position in VCF"),
            }
        }

        Ok(())
    }
}
