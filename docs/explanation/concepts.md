---
title: Core Concepts
---

# Core Concepts

This page explains the key concepts behind clam and how they relate to population genetics analysis.

## Callable Loci

A genomic site is considered **callable** for a sample if we have sufficient confidence that we could accurately determine the genotype at that position. In practice, this means the site has adequate sequencing depth.

Sites that are not callable include:

- Regions with no coverage (unmapped reads)
- Regions with very low coverage (unreliable genotype calls)
- Regions with extremely high coverage (often repetitive elements or CNVs)

clam tracks callable sites at two levels:

1. **Sample-level**: Is this site callable for a specific individual?
2. **Site-level**: Is this site callable for the analysis (e.g., callable in enough samples)?

## Depth Thresholds

clam uses depth thresholds to determine callability. Both minimum and maximum thresholds are important:

### Minimum Depth

Sites with very low depth have unreliable genotype calls. A heterozygous site with only 2x coverage has a 50% chance of appearing homozygous simply due to sampling.

Common minimum thresholds range from 5x to 15x depending on your confidence requirements and ploidy.

### Maximum Depth

Extremely high depth often indicates:

- Repetitive regions where reads from multiple genomic locations map to the same place
- Copy number variants (CNVs) or duplications
- Systematic mapping artifacts

These regions tend to have unreliable variant calls and are typically excluded. Maximum thresholds are often set to 2-3x the mean genome-wide depth.

### Per-Chromosome Thresholds

Some chromosomes may require different thresholds:

- Sex chromosomes in samples with XY sex determination have half the expected autosomal depth
- Organellar genomes (mitochondria, chloroplasts) often have much higher depth

clam supports per-chromosome threshold files to handle these cases.

## Sample-Level vs Site-Level Filtering

clam applies thresholds in two stages:

### Sample-Level Thresholds

First, for each sample at each position, clam checks if the depth falls within the acceptable range:

```
min_depth <= sample_depth <= max_depth
```

If yes, that site is callable for that sample.

### Site-Level Thresholds

Next, clam can apply aggregate filters across all samples:

- **Proportion callable** (`-d`): What fraction of samples must be callable? Setting `-d 0.8` requires 80% of samples to pass individual thresholds.
- **Mean depth range** (`-u`, `-U`): What is the acceptable range for mean depth across all samples at a site?

These filters help identify sites that are systematically problematic across the dataset.

## Populations

Many population genetic statistics compare diversity within and between groups:

- **π (pi)**: Nucleotide diversity *within* a population
- **d~xy~**: Absolute divergence *between* two populations  
- **F~ST~**: Relative differentiation between populations

To calculate these statistics, clam needs to know which samples belong to which population. This is specified using a population file (see [Input Formats](../reference/input-formats.md)).

When populations are defined:

- `clam loci` tracks callable sites per population (how many samples in each population are callable at each site)
- `clam stat` calculates within-population (π) and between-population (d~xy~, F~ST~) statistics

Without a population file, clam treats all samples as a single population and only calculates π.

## The clam Workflow

A typical clam analysis has two main steps:

### Step 1: Generate Callable Loci (`clam loci`)

```
Depth files (D4/GVCF) → clam loci → Callable loci (Zarr)
```

This step processes depth information and applies your thresholds to determine which sites are callable. The output is a compact Zarr array storing callable counts per population at each genomic position.

### Step 2: Calculate Statistics (`clam stat`)

```
Callable loci (Zarr) + Variants (VCF) → clam stat → Statistics (TSV)
```

This step combines callable site information with variant calls to compute accurate diversity statistics. The callable loci provide the denominator (total comparisons), while the VCF provides the numerator (differences).

### Optional: Pre-collect Depth (`clam collect`)

For workflows where you want to run `loci` multiple times with different thresholds:

```
Depth files (D4/GVCF) → clam collect → Depth (Zarr) → clam loci → Callable loci (Zarr)
```

The `collect` step stores raw depth values in an efficient Zarr format. This is faster than re-reading the original depth files when testing multiple threshold configurations.

## Statistics Calculated

clam calculates the following statistics in windows across the genome:

### Nucleotide Diversity (π)

The average number of pairwise differences per site within a population:

$$\pi = \frac{\sum \text{pairwise differences}}{\sum \text{pairwise comparisons}}$$

### Absolute Divergence (d~xy~)

The average number of pairwise differences per site between two populations:

$$d_{xy} = \frac{\sum \text{between-population differences}}{\sum \text{between-population comparisons}}$$

### Fixation Index (F~ST~)

A measure of population differentiation, calculated using the Hudson estimator:

$$F_{ST} = \frac{d_{xy} - \frac{\pi_1 + \pi_2}{2}}{d_{xy}}$$

clam uses a ratio-of-averages approach, summing numerators and denominators across sites within each window before computing the final ratio.

## Runs of Homozygosity (ROH)

In populations with recent inbreeding or small effective population size, individuals may have long runs of homozygosity (ROH). When calculating diversity statistics, it can be useful to exclude samples that are within ROH regions at each site.

Non-ROH heterozygosity can serve as a proxy for the inbreeding load in a population. This is because deleterious mutations that were previously masked as heterozygotes become exposed in ROH regions, and the abundance of such mutations scales with genetic diversity ([Kardos et al. 2025](https://www.sciencedirect.com/science/article/pii/S016953472500182X)).

clam can optionally accept ROH intervals (`--roh`) and will calculate heterozygosity excluding samples in ROH regions. At each site, any sample falling within an ROH region is excluded from the heterozygosity calculation for that site.

## Heterozygosity

clam calculates heterozygosity as the proportion of heterozygous sites among callable sites:

$$\text{Heterozygosity} = \frac{\text{het\_total}}{\text{callable\_total}}$$

When callable sites are provided, clam outputs a `heterozygosity.tsv` file with heterozygosity estimates per window.

### Per-Sample Mode

When using per-sample callable masks (`--per-sample` in `clam loci`), heterozygosity is calculated for each sample individually. This provides the most accurate estimates because each sample has its own callable site mask.

For each sample in each window:

- `het_total`: Number of heterozygous genotypes for the sample
- `callable_total`: Number of sites where the sample is callable
- `heterozygosity`: `het_total / callable_total`

### Per-Population Mode

When using population-level callable counts (default `clam loci` output), heterozygosity is calculated per population by summing across all samples in the population.

For each population in each window:

- `het_total`: Sum of heterozygous genotypes across all samples in the population
- `callable_total`: Sum of callable sites across all samples in the population
- `heterozygosity`: `het_total / callable_total`

!!! note
    Per-population heterozygosity with ROH exclusion is approximate because the callable counts are aggregated at the population level. For accurate per-sample ROH exclusion, use per-sample callable masks.

### ROH-Excluded Heterozygosity

When ROH data is provided, clam also calculates heterozygosity after excluding samples in ROH regions:

- `het_not_in_roh`: Heterozygous sites where the sample is not in an ROH region
- `callable_not_in_roh`: Callable sites where the sample is not in an ROH region
- `heterozygosity_not_in_roh`: `het_not_in_roh / callable_not_in_roh`
