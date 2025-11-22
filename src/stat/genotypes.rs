use color_eyre::eyre::{bail, Context, Result};
use noodles::vcf::variant::record::samples::keys::key;
use noodles::vcf::variant::record::samples::series::Value;
use noodles::vcf::variant::record::samples::Series;
use noodles::vcf::variant::record::Samples;
use noodles::vcf::{self, Header, Record};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Genotype {
    pub ref_count: u8,
    pub alt_count: u8,
}

impl Genotype {
    #[inline]
    pub fn new(ref_count: u8, alt_count: u8) -> Self {
        Self {
            ref_count,
            alt_count,
        }
    }

    #[inline]
    pub fn is_missing(&self) -> bool {
        self.ref_count == 0 && self.alt_count == 0
    }

    #[inline]
    pub fn is_heterozygous(&self) -> bool {
        self.ref_count > 0 && self.alt_count > 0
    }

    #[inline]
    pub fn is_homozygous_ref(&self) -> bool {
        self.ref_count > 0 && self.alt_count == 0
    }

    #[inline]
    pub fn is_homozygous_alt(&self) -> bool {
        self.ref_count == 0 && self.alt_count > 0
    }

    #[inline]
    pub fn ploidy(&self) -> u8 {
        self.ref_count + self.alt_count
    }
}

fn parse_genotype_value(value: Option<Value>) -> Result<Genotype> {
    match value {
        Some(Value::Genotype(gt)) => {
            let mut ref_count = 0u8;
            let mut alt_count = 0u8;

            for result in gt.iter() {
                let (allele_position, _phasing) = result.context("Failed to parse allele")?;

                match allele_position {
                    Some(0) => ref_count += 1, // Reference allele
                    Some(1) => alt_count += 1, // First alternate allele
                    Some(_) => {}              // Other alternates (ignored for biallelic)
                    None => {}                 // Missing allele
                }
            }

            Ok(Genotype::new(ref_count, alt_count))
        }
        None => Ok(Genotype::new(0, 0)), // Missing genotype
        _ => bail!("Invalid GT value type"),
    }
}

pub fn genotypes<F>(record: &vcf::Record, header: &Header) -> Result<Vec<Genotype>> {
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
        res[sample_idx] = genotype
    }

    Ok(res)
}

pub fn fold_genotypes<T, F>(
    record: &vcf::Record,
    header: &Header,
    init: T,
    mut folder: F,
) -> Result<T>
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
