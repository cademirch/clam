---
title: Input File Formats
---

# Input File Formats

This page documents the file formats accepted by clam.

## Depth Files

### D4 Files

D4 is a compact format for storing per-base depth information.

**Formats accepted:**

- Uncompressed D4 (`.d4`)
- Bgzipped D4 with index (`.d4.gz` + `.d4.gz.gzi`)

**Generating D4 files:**

```bash
# Generate D4 from BAM using mosdepth
mosdepth --d4 sample sample.bam

# Optionally compress and index for smaller file size
bgzip sample.per-base.d4
bgzip --index sample.per-base.d4.gz
```

**Sample name extraction:**

clam extracts sample names from D4 filenames. For a file named `sample1.per-base.d4.gz`, the sample name is `sample1`.

### Merged D4 Files

A merged D4 file contains depth information for multiple samples in a single file. clam automatically detects merged D4 files and extracts all sample names from the file header.

!!! warning "Merged D4 files must not be bgzipped"
    Unlike single-sample D4 files, merged D4 files should remain uncompressed (`.d4`).

### GVCF Files

GVCF (Genomic VCF) files contain per-sample depth and genotype quality information at every position.

**Requirements:**

- Must be bgzipped (`.g.vcf.gz`)
- Must be tabix indexed (`.g.vcf.gz.tbi`)

**Generating indexed GVCFs:**

```bash
bgzip sample.g.vcf
tabix -p vcf sample.g.vcf.gz
```

**Sample name extraction:**

Sample names are extracted from the filename.

---

## VCF Files

For `clam stat`, input VCF files must be bgzipped and tabix indexed.

**Requirements:**

- Must be bgzipped (`.vcf.gz`)
- Must be tabix indexed (`.vcf.gz.tbi`)
- Should contain contig lengths in the header (recommended)

**Indexing:**

```bash
bgzip variants.vcf
tabix -p vcf variants.vcf.gz
```

---

## Population File

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
- For D4 files, sample names are derived from filenames
- For VCF/GVCF files, sample names come from the file header
- Each sample should appear exactly once
- Population names can be any string (no spaces)

---

## Sample Name File

Override automatic sample name detection from filenames. Useful when filenames contain dots (e.g., "L.sample1.d4" would be auto-detected as "L" instead of "L.sample1").

**Format:** Tab-separated, two columns, no header.

| Column | Description |
|--------|-------------|
| 1 | Filename (not full path) |
| 2 | Sample name |

**Example:**

```
L.sample1.d4.gz	L.sample1
L.sample2.d4.gz	L.sample2
sample3.d4	sample3_renamed
```

**Notes:**

- Use the filename only, not the full path
- All input files must be listed when using this option
- Sample names must be unique
- Useful when files have dots in sample identifiers (e.g., species abbreviations like "L." or "D.")
- Does not apply to multisample D4 files (sample names come from internal track names)

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
| Merged D4 | `.d4` (uncompressed only) | No | `loci`, `collect` |
| GVCF | `.g.vcf.gz` | `.g.vcf.gz.tbi` | `loci`, `collect` |
| VCF | `.vcf.gz` | `.vcf.gz.tbi` | `stat` |
| Population | `.tsv` | No | `loci`, `stat` |
| Sample names | `.tsv` | No | `loci`, `collect` |
| Chromosome list | `.txt` | No | `loci`, `stat`, `collect` |
| Thresholds | `.tsv` | No | `loci` |
| ROH | `.bed.gz` | `.bed.gz.tbi` | `stat` |
| Regions | `.bed` | No | `stat` |
