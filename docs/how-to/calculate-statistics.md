---
title: Calculate Statistics
---

# Calculate Statistics

This guide covers how to use `clam stat` to calculate population genetic statistics from a VCF using callable site information.

## Prerequisites

- A VCF file (bgzipped and tabix indexed)
- Callable loci Zarr from `clam loci` (recommended)
- Optional: population file or `--samples` TSV for d~xy~ and F~ST~ (not needed if callable zarr already contains population metadata)
- Optional: ROH file for calculating non-ROH heterozygosity

## Input Requirements

### VCF File

The input VCF must be bgzipped and tabix indexed:

```bash
bgzip variants.vcf
tabix -p vcf variants.vcf.gz
```

### Callable Sites

While optional, providing callable sites from `clam loci` is strongly recommended for accurate statistics:

```bash
clam stat -o results/ -w 10000 -c callable.zarr variants.vcf.gz
```

Without callable sites, all positions in the VCF are assumed callable, which can bias diversity estimates if you have missing data.

## Basic Usage

### Fixed-Size Windows

Calculate statistics in fixed-size windows across the genome:

```bash
clam stat -o results/ -w 10000 -c callable.zarr variants.vcf.gz
```

This calculates π (and d~xy~, F~ST~ if populations are defined) in 10kb windows.

### Custom Regions

Calculate statistics for specific regions defined in a BED file:

```bash
clam stat -o results/ --regions-file genes.bed -c callable.zarr variants.vcf.gz
```

## Specifying Populations

To calculate between-population statistics (d~xy~ and F~ST~), clam needs population assignments. There are three ways to provide them:

### Automatic from Callable Zarr (Recommended)

If populations were defined during the `clam loci` step, they are stored in the callable Zarr metadata. `clam stat` reads them automatically -- no extra flags needed:

```bash
clam stat -o results/ -w 10000 -c callable.zarr variants.vcf.gz
```

This produces `dxy.tsv` and `fst.tsv` automatically when the callable zarr contains multiple populations.

### Using --samples

Provide a samples TSV to define (or override) population assignments:

```bash
clam stat -o results/ -w 10000 -c callable.zarr -s samples.tsv variants.vcf.gz
```

Only the `sample_name` and `population` columns are used; `file_path` is ignored for `clam stat`.

### Using --population-file (Deprecated)

```bash
clam stat -o results/ -w 10000 -c callable.zarr -p populations.tsv variants.vcf.gz
```

The population file maps samples to populations:

```
sample1	PopA
sample2	PopA
sample3	PopB
sample4	PopB
```

!!! warning "Population Consistency"
    When using `--samples` or `--population-file` with a callable Zarr that has `PopulationCounts` type, the number of populations must match the number of population columns in the Zarr. If they don't match, clam will produce a clear error. Use callable sites generated with `--per-sample` mode for maximum flexibility with different population definitions at stat time.

!!! note "Sample Names"
    Sample names must exactly match those in your VCF header.

### Handling Missing Samples

By default, clam errors if the VCF contains samples not in the population file. To instead warn and exclude those samples:

```bash
clam stat -o results/ -w 10000 -c callable.zarr -p populations.tsv --force-samples variants.vcf.gz
```

## Using ROH Data

If you have runs of homozygosity (ROH) calls, clam can calculate non-ROH heterozygosity by excluding samples within ROH regions at each site:

```bash
clam stat -o results/ -w 10000 -c callable.zarr -r roh.bed.gz variants.vcf.gz
```

Non-ROH heterozygosity can serve as a proxy for the inbreeding load in a population ([Kardos et al. 2025](https://www.sciencedirect.com/science/article/pii/S016953472500182X)).

The ROH file must be:

- BED format with sample name in the 4th column
- bgzipped and tabix indexed

```bash
# Prepare ROH file
sort -k1,1 -k2,2n roh.bed > roh.sorted.bed
bgzip roh.sorted.bed
tabix -p bed roh.sorted.bed.gz
```

When ROH data is provided, the `heterozygosity.tsv` output includes additional columns for ROH-excluded statistics (`het_not_in_roh`, `callable_not_in_roh`, `heterozygosity_not_in_roh`).

See [Input Formats](../reference/input-formats.md#roh-file) for file format details.

## Filtering Chromosomes

### Exclude Specific Chromosomes

```bash
clam stat -o results/ -w 10000 -c callable.zarr -x chrM,chrY variants.vcf.gz
```

### Include Only Specific Chromosomes

```bash
clam stat -o results/ -w 10000 -c callable.zarr -i chr1,chr2,chr3 variants.vcf.gz
```

## Output Files

`clam stat` produces tab-separated files in the output directory:

| File | Description | Generated |
|------|-------------|-----------|
| `pi.tsv` | Nucleotide diversity per population | Always |
| `dxy.tsv` | Absolute divergence between populations | When >1 population |
| `fst.tsv` | Fixation index between populations | When >1 population |
| `heterozygosity.tsv` | Per-sample or per-population heterozygosity | When callable sites provided |

The `heterozygosity.tsv` file contains per-sample heterozygosity when using per-sample callable masks (`--per-sample` in `clam loci`), or per-population heterozygosity when using population counts. When ROH data is provided, additional columns report heterozygosity excluding samples in ROH regions.

See [Output Formats](../reference/output-formats.md#clam-stat-output) for column descriptions.

## Performance

Use multiple threads for faster processing:

```bash
clam stat -o results/ -w 10000 -c callable.zarr -t 16 variants.vcf.gz
```

The chunk size affects parallelization and memory usage:

```bash
clam stat -o results/ -w 10000 -c callable.zarr --chunk-size 500000 variants.vcf.gz
```

!!! tip
    If you used `clam loci` with a specific chunk size, `clam stat` automatically uses the same chunk size when reading the callable Zarr.

## Complete Example

```bash
# Calculate pi, dxy, and Fst in 50kb windows
# - Use callable sites from loci
# - Define two populations
# - Calculate non-ROH heterozygosity
# - Exclude mitochondria
# - Use 16 threads
clam stat \
  -o results/ \
  -w 50000 \
  -c callable.zarr \
  -p populations.tsv \
  -r roh.bed.gz \
  -x chrM \
  -t 16 \
  variants.vcf.gz
```

## Working with Results

### Python

```python
import pandas as pd

pi = pd.read_csv("results/pi.tsv", sep="\t")
dxy = pd.read_csv("results/dxy.tsv", sep="\t")
fst = pd.read_csv("results/fst.tsv", sep="\t")

# Plot pi across chromosome 1
chr1_pi = pi[pi["chrom"] == "chr1"]
chr1_pi.plot(x="start", y="pi", kind="line")
```

### R

```r
library(tidyverse)

pi <- read_tsv("results/pi.tsv")
dxy <- read_tsv("results/dxy.tsv")
fst <- read_tsv("results/fst.tsv")

# Plot pi by population
ggplot(pi, aes(x=start, y=pi, color=population)) +
  geom_line() +
  facet_wrap(~chrom)
```
