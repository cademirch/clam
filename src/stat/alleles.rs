use anyhow::{bail, Context, Result};
use noodles::vcf::variant::record::samples::keys::key;
use noodles::vcf::variant::record::samples::series::Value;
use noodles::vcf::variant::record::samples::Series;
use noodles::vcf::{self};
use std::iter::zip;

use super::*;

enum Homozygous {
    Ref,
    Alt,
}
enum Genotype {
    Homozygous(Homozygous),
    Heterozygous,
    Missing, // All alleles missing
}

pub struct VCFData {
    pub allele_counts: Vec<[u32; 2]>,
    pub het_counts: Vec<Vec<u32>>,
}

impl VCFData {
    pub fn new(num_samples_per_population: Vec<usize>) -> Self {
        let allele_counts: Vec<[u32; 2]> = vec![[0; 2]; num_samples_per_population.len()];
        let het_counts = num_samples_per_population
            .iter()
            .map(|count| vec![0; *count])
            .collect();

        Self {
            allele_counts,
            het_counts,
        }
    }
}

/// Count alleles in a single vcf record
pub fn count_alleles(
    record: &vcf::Record,
    header: &vcf::Header,
    population_info: &PopulationMapping,
    vcfdata: &mut VCFData,
    samples_in_roh: HashSet<String>,
) -> Result<()> {
    let samples = record.samples();
    let sample_names = header.sample_names().iter();
    let gt_series = samples
        .select(key::GENOTYPE)
        .ok_or_else(|| anyhow!("Malformed variant record: {:?}", record))?;

    for (sample_name, value) in zip(sample_names, gt_series.iter(&header)) {
        let Some(Value::Genotype(genotype)) = value? else {
            bail!(
                "Invalid GT value for sample {}, in variant record: {:?}",
                sample_name,
                record
            )
        };

        let (population_idx, internal_idx) = match population_info.lookup_sample_name(sample_name) {
            Ok(indices) => indices,
            Err(_) => {
                continue;
            }
        };

        let mut ref_count = 0;
        let mut alt_count = 0;

        for result in genotype.iter() {
            let allele = result.with_context(|| {
                format!(
                    "Failed to parse alleles for sample at {}, in variant record: {:?}",
                    sample_name, record
                )
            })?;

            let (position, _) = allele;

            match position {
                Some(0) => ref_count += 1,
                Some(1) => alt_count += 1,
                None => (),
                _ => bail!(
                    "Unexpected allele value in sample {}: {:?}",
                    sample_name,
                    position
                ),
            }
        }

        let genotype = match (ref_count, alt_count) {
            (0, 0) => Genotype::Missing,
            (0, _) => Genotype::Homozygous(Homozygous::Alt),
            (_, 0) => Genotype::Homozygous(Homozygous::Ref),
            _ => Genotype::Heterozygous,
        };

        if samples_in_roh.contains(sample_name) {
            trace!(
                "ROH sample: {}: ref:{} alt: {}",
                sample_name,
                ref_count,
                alt_count
            );

            match genotype {
                Genotype::Homozygous(Homozygous::Ref) => {
                    vcfdata.allele_counts[population_idx][0] += 1
                } // Increment reference count
                Genotype::Homozygous(Homozygous::Alt) => {
                    vcfdata.allele_counts[population_idx][1] += 1
                } // Increment alternate count
                Genotype::Heterozygous => {
                    trace!(
                        "Sample {} in ROH but heterozygous (ref:{} alt:{})",
                        sample_name,
                        ref_count,
                        alt_count
                    );
                }
                Genotype::Missing => {
                    trace!(
                        "Sample {} in ROH but missing genotype (ref:{} alt:{})",
                        sample_name,
                        ref_count,
                        alt_count
                    );
                }
            }
        } else {
            trace!("refcount: {} alt_count:{}", ref_count, alt_count);
            vcfdata.allele_counts[population_idx][0] += ref_count;
            vcfdata.allele_counts[population_idx][1] += alt_count;

            if matches!(genotype, Genotype::Heterozygous) {
                trace!("het_counts: {:?}", vcfdata.het_counts);
                vcfdata.het_counts[population_idx][internal_idx] += 1; // Increment het count outside ROH
            };
        };
    }
    trace!("Counts: {:?}", vcfdata.allele_counts);
    Ok(())
}

#[cfg(test)]
mod tests {
    use std::io::Cursor;

    use anyhow::Result;

    use noodles::vcf::io::Reader;
    use noodles::vcf::Header;

    use super::*;

    fn create_population_mapping(header: &Header, pop_data: &str) -> PopulationMapping {
        let pop_data_cursor = Cursor::new(pop_data);
        PopulationMapping::from_path(pop_data_cursor, Some(header.sample_names()))
            .expect("Failed to create PopulationMapping")
    }

    #[test]
    fn test_count_alleles_single_population() -> Result<()> {
        let vcf_data: &str = concat!(
            "##fileformat=VCFv4.3\n",
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample1\tsample2\n",
            "sq0\t1\t.\tA\tG\t.\tPASS\t.\tGT\t0/0\t0/1\n"
        );

        let pop_data: &str = "sample1\tpop0\nsample2\tpop0\n";

        let vcf_cursor = Cursor::new(vcf_data);
        let mut reader = Reader::new(vcf_cursor);
        let header = reader.read_header()?;

        let population_mapping = create_population_mapping(&header, pop_data);

        for result in reader.records() {
            let record = result?;
            let mut counts = VCFData::new(population_mapping.get_sample_counts_per_population());

            count_alleles(
                &record,
                &header,
                &population_mapping,
                &mut counts,
                HashSet::new(),
            )?;

            assert_eq!(counts.allele_counts[0][0], 3); // 2 from sample1, 1 from sample2
            assert_eq!(counts.allele_counts[0][1], 1); // 1 from sample2
        }

        Ok(())
    }

    #[test]
    fn test_count_alleles_haploid_single_population() -> Result<()> {
        let vcf_data: &str = concat!(
            "##fileformat=VCFv4.3\n",
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample1\tsample2\n",
            "sq0\t1\t.\tA\tG\t.\tPASS\t.\tGT\t0\t1\n"
        );

        let pop_data: &str = "sample1\tpop0\nsample2\tpop0\n";

        let vcf_cursor = Cursor::new(vcf_data);
        let mut reader = Reader::new(vcf_cursor);
        let header = reader.read_header()?;

        let population_mapping = create_population_mapping(&header, pop_data);

        for result in reader.records() {
            let record = result?;
            trace!(
                "{:?}",
                population_mapping.get_sample_counts_per_population()
            );
            let mut counts = VCFData::new(population_mapping.get_sample_counts_per_population());

            count_alleles(
                &record,
                &header,
                &population_mapping,
                &mut counts,
                HashSet::new(),
            )?;

            assert_eq!(counts.allele_counts[0][0], 1); // 1 from sample1
            assert_eq!(counts.allele_counts[0][1], 1); // 1 from sample2
        }

        Ok(())
    }

    #[test]
    fn test_count_alleles_two_populations() -> Result<()> {
        let vcf_data: &str = concat!(
            "##fileformat=VCFv4.3\n",
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample1\tsample2\n",
            "sq0\t1\t.\tA\tG\t.\tPASS\t.\tGT\t0/0\t0/1\n"
        );

        let pop_data: &str = "sample1\tpop0\nsample2\tpop1\n";

        let vcf_cursor = Cursor::new(vcf_data);
        let mut reader = Reader::new(vcf_cursor);
        let header = reader.read_header()?;

        let population_mapping = create_population_mapping(&header, pop_data);

        for result in reader.records() {
            let record = result?;
            let mut counts = VCFData::new(population_mapping.get_sample_counts_per_population());
            trace!(
                "{:?}",
                population_mapping.get_sample_counts_per_population()
            );
            count_alleles(
                &record,
                &header,
                &population_mapping,
                &mut counts,
                HashSet::new(),
            )?;

            // Population 0
            assert_eq!(counts.allele_counts[0][0], 2); // 2 from sample1
            assert_eq!(counts.allele_counts[0][1], 0); // 0 from sample1
            assert_eq!(counts.het_counts[0][0], 0);

            // Population 1
            assert_eq!(counts.allele_counts[1][0], 1); // 1 from sample2
            assert_eq!(counts.allele_counts[1][1], 1); // 1 from sample2
            assert_eq!(counts.het_counts[1][0], 1);
        }

        Ok(())
    }

    #[test]
    fn test_count_alleles_with_roh() -> Result<()> {
        let vcf_data: &str = concat!(
            "##fileformat=VCFv4.3\n",
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample1\tsample2\n",
            "sq0\t1\t.\tA\tG\t.\tPASS\t.\tGT\t0/0\t0/1\n",
            "sq0\t2\t.\tA\tG\t.\tPASS\t.\tGT\t1/1\t0/1\n",
            "sq0\t3\t.\tA\tG\t.\tPASS\t.\tGT\t0/1\t0/1\n"
        );

        let pop_data: &str = "sample1\tpop0\nsample2\tpop0\n";

        let roh_samples: HashSet<String> = vec!["sample1".to_string()].into_iter().collect(); // Only sample1 is in ROH

        let vcf_cursor = Cursor::new(vcf_data);
        let mut reader = Reader::new(vcf_cursor);
        let header = reader.read_header()?;

        let population_mapping = create_population_mapping(&header, pop_data);

        for result in reader.records() {
            let record = result?;
            trace!(
                "{:?}",
                population_mapping.get_sample_counts_per_population()
            );
            let mut counts = VCFData::new(population_mapping.get_sample_counts_per_population());

            count_alleles(
                &record,
                &header,
                &population_mapping,
                &mut counts,
                roh_samples.clone(),
            )?;

            // Assertions depend on the record and ROH status
            match record.variant_start().unwrap().unwrap().get() {
                1 => {
                    assert_eq!(counts.allele_counts[0][0], 2); // sample1 contributes 2 ref alleles
                    assert_eq!(counts.allele_counts[0][1], 1); // sample2 contributes 1 alt allele
                }
                2 => {
                    assert_eq!(counts.allele_counts[0][0], 1); // sample2 contributes 1 ref allele
                    assert_eq!(counts.allele_counts[0][1], 2); // sample1 contributes 2 alt alleles
                }

                3 => {
                    // sample1 should be zero het counts bc in roh
                    assert_eq!(counts.het_counts[0][0], 0);
                }
                _ => panic!("Unexpected position in VCF"),
            }
        }

        Ok(())
    }
}
