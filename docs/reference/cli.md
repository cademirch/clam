---
title: CLI Reference
---

# CLI Reference

This page documents all clam commands and their options.

## clam

```
Population genetics analysis toolkit

Usage: clam <COMMAND>

Commands:
  loci     Calculate callable sites from depth statistics
  stat     Calculate population genetic statistics from VCF
  collect  Collect depth from multiple files into a Zarr store
  help     Print this message or the help of the given subcommand(s)

Options:
  -h, --help     Print help
  -V, --version  Print version
```

---

## clam loci

Calculate callable sites from depth statistics.

```
Usage: clam loci [OPTIONS] --output <OUTPUT> [INPUT]...
```

### Arguments

`[INPUT]...`
:   Input depth files. Accepts D4 files (bgzipped and indexed), GVCF files (bgzipped and tabix indexed), or a Zarr store from `clam collect`. Not required if using `--samples`. A single Zarr store can also be provided alongside `--samples` to use stored depth data with custom population assignments.

### Required Options

`-o, --output <OUTPUT>`
:   Output path for callable sites Zarr array.

### Depth Threshold Options

`-m, --min-depth <MIN_DEPTH>`
:   Minimum depth to consider a site callable for each individual. [default: 0]

`-M, --max-depth <MAX_DEPTH>`
:   Maximum depth to consider a site callable for each individual. [default: inf]

`-d, --min-proportion <MIN_PROPORTION>`
:   Proportion of samples that must pass thresholds at a site to consider it callable. Value between 0.0 and 1.0. [default: 0]

`--min-mean-depth <MEAN_DEPTH_MIN>`
:   Minimum mean depth across all samples required at a site. [default: 0]

`--max-mean-depth <MEAN_DEPTH_MAX>`
:   Maximum mean depth across all samples allowed at a site. [default: inf]

`--thresholds-file <THRESHOLD_FILE>`
:   Custom thresholds per chromosome. Tab-separated file with columns: contig, min_depth, max_depth. See [Input Formats](input-formats.md#per-chromosome-thresholds-file).

### Output Options

`--per-sample`
:   Output per-sample boolean masks instead of per-population counts. Useful when you need sample-level callability information.

### Input Options

`--min-gq <MIN_GQ>`
:   Minimum genotype quality (GQ) to count depth. Only applies to GVCF input.

### Chromosome Filtering

`-x, --exclude <EXCLUDE>...`
:   Comma-separated list of chromosomes to exclude. Example: `-x chrM,chrY`

`--exclude-file <EXCLUDE_FILE>`
:   Path to file with chromosomes to exclude, one per line.

`-i, --include <INCLUDE>...`
:   Comma-separated list of chromosomes to include (restrict analysis to). Example: `-i chr1,chr2,chr3`

`--include-file <INCLUDE_FILE>`
:   Path to file with chromosomes to include, one per line.

### Sample Input Options

`-s, --samples <SAMPLES>`
:   Path to samples TSV file specifying sample names, file paths, and optionally populations. See [Input Formats](input-formats.md#samples-file). When provided, positional input files are not allowed -- except when the sole positional argument is a Zarr store, in which case `--samples` provides population assignments only (the `file_path` column is not required).

`-p, --population-file <POPULATION_FILE>`
:   **(Deprecated)** Path to file that defines populations. Use `--samples` instead. See [Input Formats](input-formats.md#population-file).

### Performance Options

`-t, --threads <THREADS>`
:   Number of threads to use for parallel processing. [default: 1]

`--chunk-size <CHUNK_SIZE>`
:   Chunk size for processing in base pairs. [default: 1000000]

### Examples

```bash
# Basic usage with D4 files
clam loci -o callable.zarr -m 10 -M 100 sample1.d4.gz sample2.d4.gz

# Using samples file (recommended for explicit sample naming)
clam loci -o callable.zarr -m 10 --samples samples.tsv

# With population file and 80% callability threshold (deprecated)
clam loci -o callable.zarr -m 10 -d 0.8 -p populations.tsv *.d4.gz

# Using GVCF input with GQ filter
clam loci -o callable.zarr -m 10 --min-gq 20 sample1.g.vcf.gz sample2.g.vcf.gz

# Exclude sex chromosomes
clam loci -o callable.zarr -m 10 -x chrX,chrY *.d4.gz

# Using Zarr input from clam collect
clam loci -o callable.zarr -m 10 -M 100 depths.zarr

# Using Zarr input — populations are auto-read from zarr metadata
clam loci -o callable.zarr -m 10 depths.zarr

# Using Zarr input with population override via --samples
clam loci -o callable.zarr -m 10 depths.zarr --samples samples.tsv
```

---

## clam stat

Calculate population genetic statistics from a VCF using callable site information.

```
Usage: clam stat [OPTIONS] --outdir <OUTDIR> <VCF>
```

### Arguments

`<VCF>`
:   Path to input VCF file. Must be bgzipped and tabix indexed.

### Required Options

`-o, --outdir <OUTDIR>`
:   Output directory for statistics files. Will be created if it doesn't exist.

### Window Options

One of these options is required:

`-w, --window-size <WINDOW_SIZE>`
:   Window size in base pairs for calculating statistics.

`--regions-file <REGIONS_FILE>`
:   BED file specifying regions to calculate statistics for. Use this for non-overlapping custom windows or specific genomic features.

### Callable Sites Options

`-c, --callable <CALLABLE>`
:   Path to callable sites Zarr array from `clam loci`. If not provided, all sites in the VCF are assumed callable.

### ROH Options

`-r, --roh <ROH>`
:   Path to ROH (runs of homozygosity) regions BED file. Must be bgzipped and tabix indexed. Sample name should be in the 4th column. When provided, the `heterozygosity.tsv` output includes additional columns for heterozygosity calculated after excluding samples in ROH regions. See [Input Formats](input-formats.md#roh-file).

### Chromosome Filtering

`-x, --exclude <EXCLUDE>...`
:   Comma-separated list of chromosomes to exclude.

`--exclude-file <EXCLUDE_FILE>`
:   Path to file with chromosomes to exclude, one per line.

`-i, --include <INCLUDE>...`
:   Comma-separated list of chromosomes to include.

`--include-file <INCLUDE_FILE>`
:   Path to file with chromosomes to include, one per line.

### Population Options

`-s, --samples <SAMPLES>`
:   Path to samples TSV file specifying population assignments. Only `sample_name` and `population` columns are used; `file_path` is ignored. See [Input Formats](input-formats.md#samples-file).

`-p, --population-file <POPULATION_FILE>`
:   **(Deprecated)** Path to file that defines populations. Use `--samples` instead. See [Input Formats](input-formats.md#population-file).

`--force-samples`
:   Only warn about samples in VCF that are not in population file, and exclude them from analysis. Without this flag, missing samples cause an error.

!!! tip "Automatic population detection"
    When no `--samples` or `--population-file` is provided, `clam stat` automatically reads population definitions from the callable Zarr metadata (if populations were stored by `clam loci`). This means you don't need to re-specify populations if they were already defined during the `loci` step.

### Performance Options

`-t, --threads <THREADS>`
:   Number of threads to use for parallel processing. [default: 1]

`--chunk-size <CHUNK_SIZE>`
:   Chunk size for parallel processing in base pairs. [default: 1000000]

### Examples

```bash
# Basic usage with 10kb windows
clam stat -o results/ -w 10000 -c callable.zarr variants.vcf.gz

# Populations auto-read from callable zarr (no -s or -p needed)
clam stat -o results/ -w 10000 -c callable.zarr variants.vcf.gz

# With population file for dxy and Fst
clam stat -o results/ -w 10000 -c callable.zarr -p populations.tsv variants.vcf.gz

# Override populations with --samples
clam stat -o results/ -w 10000 -c callable.zarr -s samples.tsv variants.vcf.gz

# Using custom regions from BED file
clam stat -o results/ --regions-file genes.bed -c callable.zarr variants.vcf.gz

# With ROH masking
clam stat -o results/ -w 10000 -c callable.zarr -r roh.bed.gz variants.vcf.gz

# Exclude mitochondria
clam stat -o results/ -w 10000 -c callable.zarr -x chrM variants.vcf.gz
```

---

## clam collect

Collect depth from multiple files into a Zarr store. Use this when you want to run `clam loci` multiple times with different threshold parameters.

```
Usage: clam collect [OPTIONS] --output <OUTPUT> [INPUT]...
```

### Arguments

`[INPUT]...`
:   Input depth files. Accepts D4 files (bgzipped and indexed) or GVCF files (bgzipped and tabix indexed). Not required if using `--samples`.

### Required Options

`-o, --output <OUTPUT>`
:   Output path for depth Zarr array.

### Input Options

`--min-gq <MIN_GQ>`
:   Minimum genotype quality (GQ) to count depth. Only applies to GVCF input.

### Sample Input Options

`-s, --samples <SAMPLES>`
:   Path to samples TSV file specifying sample names and file paths. See [Input Formats](input-formats.md#samples-file). When provided, positional input files are not allowed.

### Chromosome Filtering

`-x, --exclude <EXCLUDE>...`
:   Comma-separated list of chromosomes to exclude.

`--exclude-file <EXCLUDE_FILE>`
:   Path to file with chromosomes to exclude, one per line.

`-i, --include <INCLUDE>...`
:   Comma-separated list of chromosomes to include.

`--include-file <INCLUDE_FILE>`
:   Path to file with chromosomes to include, one per line.

### Performance Options

`-t, --threads <THREADS>`
:   Number of threads to use for parallel processing. [default: 1]

`--chunk-size <CHUNK_SIZE>`
:   Chunk size for processing in base pairs. [default: 1000000]

### Examples

```bash
# Collect depth from D4 files
clam collect -o depths.zarr -t 8 *.d4.gz

# Using samples file (recommended for explicit sample naming)
clam collect -o depths.zarr --samples samples.tsv

# Collect from GVCFs with GQ filter
clam collect -o depths.zarr --min-gq 20 *.g.vcf.gz

# Then run loci multiple times with different thresholds
clam loci -o callable_m5.zarr -m 5 depths.zarr
clam loci -o callable_m10.zarr -m 10 depths.zarr
clam loci -o callable_m15.zarr -m 15 depths.zarr
```

---

## Shared Options

These options are available across multiple commands:

| Option | Description | Commands |
|--------|-------------|----------|
| `-t, --threads` | Number of threads for parallel processing | all |
| `-s, --samples` | Samples TSV file with sample names, file paths, and populations | `loci`, `collect`, `stat` |
| `-p, --population-file` | Population definitions file (deprecated, use `--samples`) | `loci`, `stat` |
| `-x, --exclude` | Chromosomes to exclude (comma-separated) | all |
| `--exclude-file` | File with chromosomes to exclude | all |
| `-i, --include` | Chromosomes to include (comma-separated) | all |
| `--include-file` | File with chromosomes to include | all |
| `--chunk-size` | Processing chunk size in bp | all |
