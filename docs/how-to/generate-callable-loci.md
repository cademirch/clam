---
title: Generate Callable Loci
---

# Generate Callable Loci

This guide covers how to use `clam loci` to identify genomic regions with sufficient sequencing depth to be considered callable.

## Prerequisites

- Depth files in one of the supported formats:
    - Per-sample D4 files (uncompressed or bgzipped with index)
    - Merged D4 files (uncompressed only)
    - GVCF files (bgzipped and tabix indexed)
    - A Zarr store from `clam collect`
- Optional: population file if analyzing multiple populations

## Input Formats

### D4 Files

D4 files contain per-base depth estimates. You can generate them from BAM files using [mosdepth](https://github.com/brentp/mosdepth):

```bash
# Generate D4 from BAM
mosdepth --d4 sample sample.bam

# Compress and index (required by clam)
bgzip sample.per-base.d4
bgzip --index sample.per-base.d4.gz
```

### GVCF Files

GVCF files must be bgzipped and tabix indexed:

```bash
bgzip sample.g.vcf
tabix -p vcf sample.g.vcf.gz
```

## Basic Usage

```bash
clam loci -o callable.zarr -m 10 sample1.d4.gz sample2.d4.gz sample3.d4.gz
```

This identifies sites where each sample has at least 10x depth.

## Setting Depth Thresholds

### Per-Sample Thresholds

Control which sites are callable at the individual sample level:

`-m, --min-depth`
:   Minimum depth required. Sites below this are not callable.

`-M, --max-depth`
:   Maximum depth allowed. Sites above this are not callable (often repetitive regions).

```bash
# Require 10-100x depth per sample
clam loci -o callable.zarr -m 10 -M 100 *.d4.gz
```

### Site-Level Thresholds

Apply filters across all samples at each site:

`-d, --min-proportion`
:   Fraction of samples that must pass thresholds (0.0-1.0).

`--min-mean-depth`
:   Minimum mean depth across samples.

`--max-mean-depth`
:   Maximum mean depth across samples.

```bash
# Require 80% of samples to be callable at each site
clam loci -o callable.zarr -m 10 -d 0.8 *.d4.gz
```

### Per-Chromosome Thresholds

For sex chromosomes or organellar genomes, use a thresholds file:

```bash
# thresholds.tsv
chr1	10	100
chr2	10	100
chrX	5	50
chrY	5	50
chrM	100	10000
```

```bash
clam loci -o callable.zarr --thresholds-file thresholds.tsv *.d4.gz
```

See [Input Formats](../reference/input-formats.md#per-chromosome-thresholds-file) for file format details.

## Specifying Populations

To calculate per-population callable counts (required for d~xy~ and F~ST~ downstream):

```bash
clam loci -o callable.zarr -m 10 -p populations.tsv *.d4.gz
```

The population file maps samples to populations:

```
sample1	PopA
sample2	PopA
sample3	PopB
sample4	PopB
```

!!! note "Sample Names"
    Sample names must exactly match the identifiers in your input files. For D4 files, this is typically the filename prefix (e.g., `sample1` for `sample1.d4.gz`).

See [Input Formats](../reference/input-formats.md#population-file) for details.

## Filtering Chromosomes

### Exclude Specific Chromosomes

```bash
# Inline
clam loci -o callable.zarr -m 10 -x chrM,chrY *.d4.gz

# From file
clam loci -o callable.zarr -m 10 --exclude-file exclude.txt *.d4.gz
```

### Include Only Specific Chromosomes

```bash
# Inline
clam loci -o callable.zarr -m 10 -i chr1,chr2,chr3 *.d4.gz

# From file
clam loci -o callable.zarr -m 10 --include-file autosomes.txt *.d4.gz
```

## Output Options

### Per-Sample Masks

By default, `clam loci` outputs callable counts per population. To output per-sample boolean masks instead:

```bash
clam loci -o callable.zarr -m 10 --per-sample *.d4.gz
```

This is useful when you need sample-level callability information for other analyses.

## Using GVCF Input

For GVCF files, you can filter by genotype quality (GQ):

```bash
clam loci -o callable.zarr -m 10 --min-gq 20 *.g.vcf.gz
```

## Performance

Use multiple threads for faster processing:

```bash
clam loci -o callable.zarr -m 10 -t 16 *.d4.gz
```

Adjust chunk size if needed (default 1Mb):

```bash
clam loci -o callable.zarr -m 10 --chunk-size 500000 *.d4.gz
```

## Complete Example

```bash
# Process 50 samples with populations, excluding sex chromosomes
clam loci \
  -o callable.zarr \
  -m 10 \
  -M 100 \
  -d 0.8 \
  -p populations.tsv \
  -x chrX,chrY,chrM \
  -t 16 \
  *.d4.gz
```

## Next Steps

Use the output with `clam stat` to calculate population genetic statistics:

```bash
clam stat -o results/ -w 10000 -c callable.zarr variants.vcf.gz
```

See [Calculate Statistics](calculate-statistics.md) for details.
