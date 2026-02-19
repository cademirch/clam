---
title: Input File Formats
---

# Input File Formats

This page documents the file formats accepted by clam.

## Depth Files

### D4 Files

[D4](https://github.com/38/d4-format) is a compact format for storing per-base depth information.

**Formats accepted:**

- Uncompressed D4 (`.d4`)
- Bgzipped D4 with index (`.d4.gz` + `.d4.gz.gzi`)

**Generating D4 files:**

```bash
# Generate D4 from BAM using mosdepth
mosdepth --d4 sample sample.bam

# Optionally bgzip and index to save disk space
bgzip --index sample.per-base.d4
```

**Sample name extraction:**

When using positional arguments, clam extracts sample names from D4 filenames. For a file named `sample1.per-base.d4.gz`, the sample name is `sample1` (the part before the first `.`).

To use explicit sample names instead, use the `--samples` option with a [samples file](#samples-file).

### GVCF Files

[GVCF (Genomic VCF)](https://github.com/broadinstitute/gatk-docs/blob/master/gatk3-faqs/What_is_a_GVCF_and_how_is_it_different_from_a_%27regular%27_VCF%3F.md) files contain per-sample depth and genotype quality information at every position.

**Requirements:**

- Must be bgzipped (`.g.vcf.gz`)
- Must be tabix indexed (`.g.vcf.gz.tbi`)

**Generating indexed GVCFs:**

```bash
bgzip sample.g.vcf
tabix -p vcf sample.g.vcf.gz
```

**Sample name extraction:**

When using positional arguments, sample names are extracted from the filename. For a file named `sample1.g.vcf.gz`, the sample name is `sample1` (the part before the first `.`).

To use explicit sample names instead, use the `--samples` option with a [samples file](#samples-file).

---

## VCF Files

For `clam stat`, input VCF files must be bgzipped and tabix indexed.

**Requirements:**

- Must be bgzipped (`.vcf.gz`)
- Must be tabix indexed (`.vcf.gz.tbi`)

---

## Samples File

Specifies sample names, input file paths, and optionally population assignments. This is the recommended way to provide input to `clam loci` and `clam collect` as it allows explicit control over sample naming.

**Format:** Tab-separated with header. The `sample_name` and `file_path` columns are required; `population` is optional.

| Column        | Required | Description                     |
| ------------- | -------- | ------------------------------- |
| `sample_name` | Yes      | Unique sample identifier        |
| `file_path`   | Conditional | Path to depth file (D4 or GVCF). Required for `collect` and `loci` with raw depth input. Not required when used with `loci` + zarr input (population override) or with `stat`. |
| `population`  | No       | Population assignment           |

**Example with populations:**

```
sample_name	file_path	population
sample1	/path/to/sample1.d4.gz	PopA
sample2	/path/to/sample2.d4.gz	PopA
sample3	/path/to/sample3.d4.gz	PopB
sample4	/path/to/sample4.d4.gz	PopB
```

**Example without populations:**

```
sample_name	file_path
sample1	/path/to/abc.sample1.d4.gz
sample2	/path/to/abc.sample2.d4.gz
sample3	/path/to/abc.sample3.d4.gz
```

**Example for population override (no file_path):**

```
sample_name	population
sample1	PopA
sample2	PopA
sample3	PopB
sample4	PopB
```

This format is used when providing `--samples` to:

- `clam loci` with a Zarr positional argument (population override only)
- `clam stat` (population assignment only)

**Notes:**

- Column order doesn't matter (detected from header)
- Sample names must be unique
- If the `population` column is present, all rows must have values
- If the `population` column is absent, all samples are assigned to a "default" population
- File paths can be absolute or relative to the current working directory
- When `file_path` is not required (e.g., `stat` or `loci` with zarr input), the column can be omitted entirely
- Use this format when filenames don't match desired sample names (e.g., files named `abc.sample1.d4.gz` but you want sample name `sample1`)

**Usage:**

```bash
clam loci -o callable.zarr -m 10 --samples samples.tsv
clam loci -o callable.zarr -m 10 depths.zarr --samples samples.tsv  # population override
clam collect -o depths.zarr --samples samples.tsv
clam stat -o results/ -w 10000 -c callable.zarr --samples samples.tsv
```

---

## Population File

!!! warning "Deprecated"
    The `--population-file` option is deprecated. Use `--samples` with a `population` column instead.

Defines which samples belong to which population. Used by both `clam loci` and `clam stat`.

**Format:** Tab-separated, two columns, no header.

| Column | Description |
|--------|-------------|
| 1 | Sample name |
| 2 | Population name |

**Example:**

```
sample1	PopA
sample2	PopA
sample3	PopA
sample4	PopB
sample5	PopB
sample6	PopC
```

**Notes:**

- Sample names must exactly match those in your input files
- For D4 files, sample names are derived from filenames (the part before the first `.`)
- For VCF/GVCF files, sample names come from the file header
- Each sample should appear exactly once
- Population names can be any string (no spaces)

---

## Chromosome Include/Exclude Files

Specify chromosomes to include or exclude from analysis.

**Format:** One chromosome name per line, no header.

**Example (`exclude_chroms.txt`):**

```
chrM
chrY
chrUn_gl000220
```

**Example (`include_chroms.txt`):**

```
chr1
chr2
chr3
chr4
chr5
```

**Usage:**

```bash
# Exclude specific chromosomes
clam loci --exclude-file exclude_chroms.txt ...

# Only analyze specific chromosomes
clam loci --include-file include_chroms.txt ...
```

---

## Per-Chromosome Thresholds File

Specify different depth thresholds for different chromosomes. Useful for sex chromosomes or organellar genomes.

**Format:** Tab-separated, three columns, no header.

| Column | Description |
|--------|-------------|
| 1 | Chromosome name |
| 2 | Minimum depth |
| 3 | Maximum depth |

**Example:**

```
chr1	10	100
chr2	10	100
chrX	5	50
chrY	5	50
chrM	100	10000
```

**Notes:**

- Chromosomes not listed in the file will use the default thresholds from `-m` and `-M` options
- This allows setting lower thresholds for hemizygous chromosomes (X, Y in XY individuals)
- Mitochondrial/chloroplast genomes often need much higher thresholds

---

## ROH File

Specifies runs of homozygosity (ROH) regions per sample. Used by `clam stat` to exclude samples within ROH regions when calculating diversity, enabling estimation of non-ROH heterozygosity (π~non-ROH~).

**Format:** BED format (tab-separated) with sample name in the 4th column.

| Column | Description |
|--------|-------------|
| 1 | Chromosome |
| 2 | Start position (0-based) |
| 3 | End position |
| 4 | Sample name |

**Requirements:**

- Must be bgzipped (`.bed.gz`)
- Must be tabix indexed (`.bed.gz.tbi`)

**Example (`roh.bed`):**

```
chr1	1000000	2000000	sample1
chr1	5000000	5500000	sample1
chr1	1500000	2500000	sample2
chr2	3000000	4000000	sample1
```

**Preparing the file:**

```bash
# Sort by chromosome and position
sort -k1,1 -k2,2n roh.bed > roh.sorted.bed

# Compress and index
bgzip roh.sorted.bed
tabix -p bed roh.sorted.bed.gz
```

**Notes:**

- Sample names must match those in your VCF
- ROH regions can overlap between samples
- Coordinates are 0-based, half-open (standard BED format)

---

## Regions File

Specifies custom regions for calculating statistics in `clam stat`. Use this instead of `--window-size` for non-uniform windows or specific genomic features.

**Format:** Standard BED format (tab-separated).

| Column | Description |
|--------|-------------|
| 1 | Chromosome |
| 2 | Start position (0-based) |
| 3 | End position |

**Example (`genes.bed`):**

```
chr1	11869	14409
chr1	14404	29570
chr1	17369	17436
chr2	38814	46588
```

**Notes:**

- Coordinates are 0-based, half-open (standard BED format)
- Regions can be of variable size
- Regions should not overlap (statistics are calculated independently per region)
- Additional columns (name, score, etc.) are ignored

---

## Summary Table

| File Type | Extension | Index Required | Used By |
|-----------|-----------|----------------|---------|
| D4 depth | `.d4` or `.d4.gz` | `.d4.gz.gzi` (if compressed) | `loci`, `collect` |
| GVCF | `.g.vcf.gz` | `.g.vcf.gz.tbi` | `loci`, `collect` |
| VCF | `.vcf.gz` | `.vcf.gz.tbi` | `stat` |
| Samples | `.tsv` | No | `loci`, `collect`, `stat` |
| Population | `.tsv` | No | `loci`, `stat` (deprecated for `loci`) |
| Chromosome list | `.txt` | No | `loci`, `stat`, `collect` |
| Thresholds | `.tsv` | No | `loci` |
| ROH | `.bed.gz` | `.bed.gz.tbi` | `stat` |
| Regions | `.bed` | No | `stat` |
