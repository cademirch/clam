use color_eyre::eyre::{bail, Context, Result};
use ndarray::{Array1, Array2};
use noodles::vcf::variant::record::samples::keys::key;
use noodles::vcf::variant::record::samples::series::Value;
use noodles::vcf::variant::record::samples::Series;
use noodles::vcf::variant::record::Samples;
use noodles::vcf::variant::RecordBuf;
use noodles::vcf::{Header, Record};
use std::ops::{Add, AddAssign};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[repr(u8)]
pub enum Ploidy {
    Haploid = 1,
    Diploid = 2,
}

impl Ploidy {
    #[inline]
    pub fn as_u8(self) -> u8 {
        self as u8
    }

    #[inline]
    pub fn as_usize(self) -> usize {
        self as usize
    }

    /// Infer ploidy from the first sample in a VCF record
    pub fn from_record_buf(record: &RecordBuf, header: &Header) -> Result<Self> {
        let samples = record.samples();
        let gt_series = samples
            .select(key::GENOTYPE)
            .ok_or_else(|| color_eyre::eyre::eyre!("No GT field in variant record"))?;

        let first_sample = gt_series
            .iter(header)
            .next()
            .ok_or_else(|| color_eyre::eyre::eyre!("No samples in record"))??;

        match first_sample {
            Some(Value::Genotype(genotype)) => {
                let mut ploidy = 0u8;
                for result in genotype.iter() {
                    let (_allele_position, _phasing) = result.context("Failed to parse allele")?;
                    ploidy += 1;
                }

                match ploidy {
                    1 => Ok(Ploidy::Haploid),
                    2 => Ok(Ploidy::Diploid),
                    _ => bail!("Unsupported ploidy: {}", ploidy),
                }
            }
            None => {
                // Missing genotype in haploid is represented as None
                Ok(Ploidy::Haploid)
            }
            _ => bail!("Invalid GT value type"),
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Genotype {
    pub refs: u8,
    pub alts: u8,
}

impl Genotype {
    #[inline]
    pub fn is_missing(self) -> bool {
        self.refs == 0 && self.alts == 0
    }

    #[inline]
    pub fn is_het(self) -> bool {
        self.refs > 0 && self.alts > 0
    }
}

fn parse_genotype_value(value: Option<Value>) -> Result<Genotype> {
    match value {
        Some(Value::Genotype(gt)) => {
            let mut refs = 0u8;
            let mut alts = 0u8;

            for result in gt.iter() {
                let (allele_position, _phasing) = result.context("Failed to parse allele")?;

                match allele_position {
                    Some(0) => refs += 1,
                    Some(1) => alts += 1,
                    Some(_) => {} // ignore other alleles
                    None => {}    // ignore missing alleles
                }
            }

            Ok(Genotype { refs, alts })
        }
        None => Ok(Genotype { refs: 0, alts: 0 }),
        _ => bail!("Invalid GT value type"),
    }
}

pub fn genotypes(record: &Record, header: &Header) -> Result<Vec<Genotype>> {
    let samples = record.samples();
    let mut res: Vec<Genotype> = Vec::with_capacity(samples.len());
    let gt_series = samples
        .select(key::GENOTYPE)
        .ok_or_else(|| color_eyre::eyre::eyre!("No GT field in variant record"))?;

    for (sample_idx, result) in gt_series.iter(header).enumerate() {
        let value = result.map_err(|e| {
            color_eyre::eyre::eyre!("Failed to read GT value at sample {}: {}", sample_idx, e)
        })?;
        let genotype = parse_genotype_value(value)?;
        res.push(genotype);
    }

    Ok(res)
}

pub fn fold_genotypes<T, F>(record: Record, header: &Header, init: T, mut folder: F) -> Result<T>
where
    F: FnMut(T, usize, Genotype) -> T,
{
    let samples = record.samples();
    let gt_series = samples
        .select(key::GENOTYPE)
        .ok_or_else(|| color_eyre::eyre::eyre!("No GT field in variant record"))?;

    let mut acc = init;
    for (sample_idx, result) in gt_series.iter(header).enumerate() {
        let value = result.map_err(|e| {
            color_eyre::eyre::eyre!("Failed to read GT value at sample {}: {}", sample_idx, e)
        })?;
        let genotype = parse_genotype_value(value)?;
        acc = folder(acc, sample_idx, genotype);
    }

    Ok(acc)
}

/// Represents counts of alleles (ref/alt) and missing genotypes computed from genotypes at a single variant position
#[derive(Debug, Clone)]
pub struct AlleleCounts {
    /// Reference allele count per population (ALL samples)
    pub refs: Array1<usize>,

    /// Alternate allele count per population (ALL samples)
    pub alts: Array1<usize>,

    /// Missing genotype count per population (ALL samples)
    pub missing: Array1<usize>,

    /// Heterozygous sites per sample
    pub hets: Array1<bool>,

    /// Counts excluding ROH samples (None if no ROH data)
    pub non_roh: Option<NonRohCounts>,
}

/// Allele counts excluding samples in runs of homozygosity
#[derive(Debug, Clone)]
pub struct NonRohCounts {
    /// Reference allele count per population (NON-ROH samples only)
    pub refs: Array1<usize>,

    /// Alternate allele count per population (NON-ROH samples only)
    pub alts: Array1<usize>,

    /// Heterozygous sites per sample (NON-ROH samples only)
    pub hets: Array1<bool>,
}

impl AlleleCounts {
    pub fn from_genotypes(
        genotypes: &[Genotype],
        ploidy: Ploidy,
        samples_in_roh: Option<&[usize]>,
        pop_membership: &Array2<usize>,
    ) -> Result<Self> {
        let num_samples = genotypes.len();
        let ploidy_usize = ploidy.as_usize();

        let ref_alleles = Array1::from_shape_fn(num_samples, |i| genotypes[i].refs as usize);
        let alt_alleles = Array1::from_shape_fn(num_samples, |i| genotypes[i].alts as usize);
        let missing = Array1::from_shape_fn(num_samples, |i| {
            genotypes[i].is_missing() as usize * ploidy_usize
        });
        let hets = Array1::from_shape_fn(num_samples, |i| genotypes[i].is_het());

        let refs = ref_alleles.dot(pop_membership);
        let alts = alt_alleles.dot(pop_membership);
        let missing = missing.dot(pop_membership);

        let non_roh = samples_in_roh.map(|roh_samples| {
            let mut ref_alleles_non_roh = ref_alleles.clone();
            let mut alt_alleles_non_roh = alt_alleles.clone();
            let mut hets_non_roh = hets.clone();

            // Zero out ROH samples
            for &sample_idx in roh_samples {
                ref_alleles_non_roh[sample_idx] = 0;
                alt_alleles_non_roh[sample_idx] = 0;
                hets_non_roh[sample_idx] = false;
            }

            // Zero out missing genotypes
            for (i, gt) in genotypes.iter().enumerate() {
                if gt.is_missing() {
                    ref_alleles_non_roh[i] = 0;
                    alt_alleles_non_roh[i] = 0;
                    hets_non_roh[i] = false;
                }
            }

            NonRohCounts {
                refs: ref_alleles_non_roh.dot(pop_membership),
                alts: alt_alleles_non_roh.dot(pop_membership),
                hets: hets_non_roh,
            }
        });

        Ok(Self {
            refs,
            alts,
            missing,
            hets,
            non_roh,
        })
    }
}

impl Add for AlleleCounts {
    type Output = Self;

    fn add(self, other: Self) -> Self {
        Self {
            refs: &self.refs + &other.refs,
            alts: &self.alts + &other.alts,
            missing: &self.missing + &other.missing,
            hets: self.hets | other.hets, // logical OR - a sample is het if it's het at any position
            non_roh: match (self.non_roh, other.non_roh) {
                (Some(a), Some(b)) => Some(a + b),
                _ => None,
            },
        }
    }
}

impl Add for NonRohCounts {
    type Output = Self;

    fn add(self, other: Self) -> Self {
        Self {
            refs: &self.refs + &other.refs,
            alts: &self.alts + &other.alts,
            hets: self.hets | other.hets,
        }
    }
}

impl AddAssign for AlleleCounts {
    fn add_assign(&mut self, other: Self) {
        self.refs += &other.refs;
        self.alts += &other.alts;
        self.missing += &other.missing;
        self.hets = &self.hets | &other.hets;

        if let (Some(self_non_roh), Some(other_non_roh)) = (&mut self.non_roh, other.non_roh) {
            *self_non_roh += other_non_roh;
        }
    }
}

impl AddAssign for NonRohCounts {
    fn add_assign(&mut self, other: Self) {
        self.refs += &other.refs;
        self.alts += &other.alts;
        self.hets = &self.hets | &other.hets;
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core::population::PopulationMap;
    use indexmap::IndexMap;
    use ndarray::array;
    use noodles::vcf;

    fn create_test_header() -> Header {
        let mut header = Header::default();
        header.sample_names_mut().insert("sample1".to_string());
        header.sample_names_mut().insert("sample2".to_string());
        header.sample_names_mut().insert("sample3".to_string());
        header
    }

    fn create_test_population_map() -> PopulationMap {
        let mut pop_data = IndexMap::new();
        pop_data.insert(
            "pop0".to_string(),
            vec!["sample1".to_string(), "sample3".to_string()],
        );
        pop_data.insert("pop1".to_string(), vec!["sample2".to_string()]);

        PopulationMap::from_populations(pop_data).expect("Valid population map")
    }

    fn create_test_samples() -> Vec<String> {
        vec![
            "sample1".to_string(),
            "sample2".to_string(),
            "sample3".to_string(),
        ]
    }

    fn parse_vcf_line(line: &str, header: &Header) -> Result<Record> {
        let mut reader = vcf::io::Reader::new(line.as_bytes());
        let mut record = vcf::Record::default();
        reader.read_record(&mut record)?;
        Ok(record)
    }

    #[test]
    fn test_ploidy_as_u8() {
        assert_eq!(Ploidy::Haploid.as_u8(), 1);
        assert_eq!(Ploidy::Diploid.as_u8(), 2);
    }

    #[test]
    fn test_ploidy_from_record_diploid() -> Result<()> {
        let header = create_test_header();
        let line = "chr1\t100\t.\tA\tT\t.\tPASS\t.\tGT\t0/1\t1/1\t0/0\n";
        let record = parse_vcf_line(line, &header)?;
        let record_buf = super::RecordBuf::try_from_variant_record(&header, &record)?;
        let ploidy = Ploidy::from_record_buf(&record_buf, &header)?;
        assert_eq!(ploidy, Ploidy::Diploid);
        Ok(())
    }

    #[test]
    fn test_ploidy_from_record_haploid() -> Result<()> {
        let header = create_test_header();
        let line = "chr1\t100\t.\tA\tT\t.\tPASS\t.\tGT\t0\t1\t0\n";
        let record = parse_vcf_line(line, &header)?;

        let record_buf = super::RecordBuf::try_from_variant_record(&header, &record)?;
        let ploidy = Ploidy::from_record_buf(&record_buf, &header)?;
        assert_eq!(ploidy, Ploidy::Haploid);
        Ok(())
    }

    #[test]
    fn test_ploidy_from_record_missing() -> Result<()> {
        let header = create_test_header();
        let line = "chr1\t100\t.\tA\tT\t.\tPASS\t.\tGT\t.\t.\t.\n";
        let record = parse_vcf_line(line, &header)?;

        let record_buf = super::RecordBuf::try_from_variant_record(&header, &record)?;
        let ploidy = Ploidy::from_record_buf(&record_buf, &header)?;
        assert_eq!(ploidy, Ploidy::Haploid);
        Ok(())
    }

    #[test]
    fn test_genotypes_diploid() -> Result<()> {
        let header = create_test_header();
        let line = "chr1\t100\t.\tA\tT\t.\tPASS\t.\tGT\t0/0\t0/1\t1/1\n";
        let record = parse_vcf_line(line, &header)?;

        let gts = genotypes(&record, &header)?;

        assert_eq!(gts.len(), 3);
        assert_eq!(gts[0], Genotype { refs: 2, alts: 0 });
        assert_eq!(gts[1], Genotype { refs: 1, alts: 1 });
        assert_eq!(gts[2], Genotype { refs: 0, alts: 2 });
        Ok(())
    }

    #[test]
    fn test_genotypes_haploid() -> Result<()> {
        let header = create_test_header();
        let line = "chr1\t100\t.\tA\tT\t.\tPASS\t.\tGT\t0\t1\t0\n";
        let record = parse_vcf_line(line, &header)?;

        let gts = genotypes(&record, &header)?;

        assert_eq!(gts.len(), 3);
        assert_eq!(gts[0], Genotype { refs: 1, alts: 0 });
        assert_eq!(gts[1], Genotype { refs: 0, alts: 1 });
        assert_eq!(gts[2], Genotype { refs: 1, alts: 0 });
        Ok(())
    }

    #[test]
    fn test_genotypes_with_missing() -> Result<()> {
        let header = create_test_header();
        let line = "chr1\t100\t.\tA\tT\t.\tPASS\t.\tGT\t0/0\t./.\t1/1\n";
        let record = parse_vcf_line(line, &header)?;

        let gts = genotypes(&record, &header)?;

        assert_eq!(gts.len(), 3);
        assert_eq!(gts[0], Genotype { refs: 2, alts: 0 });
        assert_eq!(gts[1], Genotype { refs: 0, alts: 0 });
        assert_eq!(gts[2], Genotype { refs: 0, alts: 2 });
        Ok(())
    }

    #[test]
    fn test_genotypes_phased() -> Result<()> {
        let header = create_test_header();
        let line = "chr1\t100\t.\tA\tT\t.\tPASS\t.\tGT\t0|0\t0|1\t1|1\n";
        let record = parse_vcf_line(line, &header)?;

        let gts = genotypes(&record, &header)?;

        assert_eq!(gts.len(), 3);
        assert_eq!(gts[0], Genotype { refs: 2, alts: 0 });
        assert_eq!(gts[1], Genotype { refs: 1, alts: 1 });
        assert_eq!(gts[2], Genotype { refs: 0, alts: 2 });
        Ok(())
    }

    #[test]
    fn test_genotypes_partial_missing() -> Result<()> {
        let header = create_test_header();
        let line = "chr1\t100\t.\tA\tT\t.\tPASS\t.\tGT\t./0\t./1\t0/.\n";
        let record = parse_vcf_line(line, &header)?;

        let gts = genotypes(&record, &header)?;

        assert_eq!(gts.len(), 3);
        // ./0 -> only ref allele observed
        assert_eq!(gts[0], Genotype { refs: 1, alts: 0 });
        // ./1 -> only alt allele observed
        assert_eq!(gts[1], Genotype { refs: 0, alts: 1 });
        // 0/. -> only ref allele observed
        assert_eq!(gts[2], Genotype { refs: 1, alts: 0 });
        Ok(())
    }

    #[test]
    fn test_genotype_is_het() {
        assert!(!Genotype { refs: 0, alts: 0 }.is_het()); // missing
        assert!(!Genotype { refs: 2, alts: 0 }.is_het()); // homozygous ref
        assert!(!Genotype { refs: 0, alts: 2 }.is_het()); // homozygous alt
        assert!(Genotype { refs: 1, alts: 1 }.is_het()); // heterozygous
        assert!(!Genotype { refs: 1, alts: 0 }.is_het()); // partial missing ./0
    }

    #[test]
    fn test_genotype_is_missing() {
        assert!(Genotype { refs: 0, alts: 0 }.is_missing());
        assert!(!Genotype { refs: 2, alts: 0 }.is_missing());
        assert!(!Genotype { refs: 0, alts: 2 }.is_missing());
        assert!(!Genotype { refs: 1, alts: 1 }.is_missing());
        assert!(!Genotype { refs: 1, alts: 0 }.is_missing()); // partial missing ./0 is NOT fully missing
    }

    #[test]
    fn test_allele_counts_no_roh() -> Result<()> {
        let genotypes = vec![
            Genotype { refs: 2, alts: 0 }, // HomozygousRef
            Genotype { refs: 1, alts: 1 }, // Heterozygous
            Genotype { refs: 0, alts: 2 }, // HomozygousAlt
        ];

        let pop_map = create_test_population_map();
        let samples = create_test_samples();
        let pop_membership = pop_map.membership_matrix(&samples)?;

        let counts =
            AlleleCounts::from_genotypes(&genotypes, Ploidy::Diploid, None, &pop_membership)?;

        // Pop 0: sample1 (2 ref) + sample3 (0 ref) = 2
        // Pop 1: sample2 (1 ref) = 1
        assert_eq!(counts.refs, array![2, 1]);

        // Pop 0: sample1 (0 alt) + sample3 (2 alt) = 2
        // Pop 1: sample2 (1 alt) = 1
        assert_eq!(counts.alts, array![2, 1]);

        // No missing genotypes
        assert_eq!(counts.missing, array![0, 0]);

        // Hets: only sample2
        assert_eq!(counts.hets, array![false, true, false]);

        assert!(counts.non_roh.is_none());

        Ok(())
    }

    #[test]
    fn test_allele_counts_with_missing() -> Result<()> {
        let genotypes = vec![
            Genotype { refs: 2, alts: 0 }, // HomozygousRef
            Genotype { refs: 0, alts: 0 }, // Missing
            Genotype { refs: 0, alts: 2 }, // HomozygousAlt
        ];

        let pop_map = create_test_population_map();
        let samples = create_test_samples();
        let pop_membership = pop_map.membership_matrix(&samples)?;

        let counts =
            AlleleCounts::from_genotypes(&genotypes, Ploidy::Diploid, None, &pop_membership)?;

        // Pop 0: sample1 (2 ref) + sample3 (0 ref) = 2
        // Pop 1: sample2 (missing) = 0
        assert_eq!(counts.refs, array![2, 0]);

        // Pop 0: sample1 (0 alt) + sample3 (2 alt) = 2
        // Pop 1: sample2 (missing) = 0
        assert_eq!(counts.alts, array![2, 0]);

        // Pop 0: 0 missing
        // Pop 1: sample2 is missing = 2 (diploid missing counts as 2)
        assert_eq!(counts.missing, array![0, 2]);

        Ok(())
    }

    #[test]
    fn test_allele_counts_with_missing_haploid() -> Result<()> {
        let genotypes = vec![
            Genotype { refs: 1, alts: 0 }, // HomozygousRef
            Genotype { refs: 0, alts: 0 }, // Missing
            Genotype { refs: 0, alts: 1 }, // HomozygousAlt
        ];

        let pop_map = create_test_population_map();
        let samples = create_test_samples();
        let pop_membership = pop_map.membership_matrix(&samples)?;

        let counts =
            AlleleCounts::from_genotypes(&genotypes, Ploidy::Haploid, None, &pop_membership)?;

        // Missing in haploid counts as 1
        assert_eq!(counts.missing, array![0, 1]);

        Ok(())
    }

    #[test]
    fn test_allele_counts_with_roh() -> Result<()> {
        let genotypes = vec![
            Genotype { refs: 2, alts: 0 }, // HomozygousRef
            Genotype { refs: 1, alts: 1 }, // Heterozygous
            Genotype { refs: 0, alts: 2 }, // HomozygousAlt
        ];

        // Samples 0 and 2 are in ROH
        let samples_in_roh = vec![0, 2];

        let pop_map = create_test_population_map();
        let samples = create_test_samples();
        let pop_membership = pop_map.membership_matrix(&samples)?;

        let counts = AlleleCounts::from_genotypes(
            &genotypes,
            Ploidy::Diploid,
            Some(&samples_in_roh),
            &pop_membership,
        )?;

        // All samples counts
        assert_eq!(counts.refs, array![2, 1]);
        assert_eq!(counts.alts, array![2, 1]);

        let non_roh = counts.non_roh.as_ref().unwrap();

        // NON-ROH counts should only include sample2 (in pop1)
        // sample1 and sample3 are in ROH, so excluded
        assert_eq!(non_roh.refs, array![0, 1]);
        assert_eq!(non_roh.alts, array![0, 1]);

        // Only sample2 is het and not in ROH
        assert_eq!(non_roh.hets, array![false, true, false]);

        Ok(())
    }

    #[test]
    fn test_allele_counts_non_roh_excludes_missing() -> Result<()> {
        let genotypes = vec![
            Genotype { refs: 2, alts: 0 }, // HomozygousRef
            Genotype { refs: 0, alts: 0 }, // Missing
            Genotype { refs: 0, alts: 2 }, // HomozygousAlt
        ];

        // Sample 0 is in ROH
        let samples_in_roh = vec![0];

        let pop_map = create_test_population_map();
        let samples = create_test_samples();
        let pop_membership = pop_map.membership_matrix(&samples)?;

        let counts = AlleleCounts::from_genotypes(
            &genotypes,
            Ploidy::Diploid,
            Some(&samples_in_roh),
            &pop_membership,
        )?;

        let non_roh = counts.non_roh.as_ref().unwrap();

        // NON-ROH counts should only include sample3 (sample1 in ROH, sample2 missing)
        // Pop 0: sample3 (0 ref) = 0
        // Pop 1: sample2 excluded (missing) = 0
        assert_eq!(non_roh.refs, array![0, 0]);
        assert_eq!(non_roh.alts, array![2, 0]);

        Ok(())
    }

    #[test]
    fn test_allele_counts_het_excludes_missing() -> Result<()> {
        let genotypes = vec![
            Genotype { refs: 1, alts: 1 }, // Heterozygous
            Genotype { refs: 0, alts: 0 }, // Missing
            Genotype { refs: 1, alts: 1 }, // Heterozygous
        ];

        // No samples in ROH
        let samples_in_roh: Vec<usize> = vec![];

        let pop_map = create_test_population_map();
        let samples = create_test_samples();
        let pop_membership = pop_map.membership_matrix(&samples)?;

        let counts = AlleleCounts::from_genotypes(
            &genotypes,
            Ploidy::Diploid,
            Some(&samples_in_roh),
            &pop_membership,
        )?;

        let non_roh = counts.non_roh.as_ref().unwrap();

        // Het counts should exclude missing genotypes
        assert_eq!(non_roh.hets, array![true, false, true]);

        Ok(())
    }
}
