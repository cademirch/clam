---
title: Output File Formats
---

# Output File Formats

This page documents the output files produced by clam commands.

## Zarr Format Overview

clam uses [Zarr](https://zarr.dev/) for its intermediate and output files. Zarr is a format for storing chunked, compressed N-dimensional arrays.

!!! warning "Schema Stability"
    The Zarr schema (metadata structure, array layout) is not yet stable and may change in future versions of clam. If you are building tools that read clam's Zarr outputs, be prepared to update them when upgrading clam.

### Why Zarr?

**Efficient Compression**
:   clam's Zarr outputs use Zstd compression with byte shuffling, typically achieving 5-10x compression ratios. The Zarr output from `clam collect` is often smaller than the equivalent merged D4 file.

**Chunked Storage**
:   Zarr stores data in chunks, enabling parallel processing, random access to specific genomic regions, and streaming without loading everything into memory.

**Interoperability**
:   Zarr is widely supported: Python (`zarr`, `xarray`, `dask`), R (`pizzarr`), Julia (`Zarr.jl`), and JavaScript (`zarr.js`). Load clam outputs directly into your analysis environment without format conversion.

**Self-Describing**
:   Zarr files include metadata describing the array structure, dimensions, and custom attributes. clam stores contig information, sample names, and chunk sizes in the metadata.

### Compression Details

**Numeric Data (depth, counts)**

- Codec: Blosc with Zstd compressor
- Compression level: 5
- Shuffle: Byte shuffle (improves compression of numeric data)

**Boolean Data (per-sample masks)**

- Codec: PackBits (efficient for boolean arrays)

---

## clam collect Output

### Depth Zarr Store

`clam collect` produces a Zarr store containing raw depth values for all samples.

**Structure:**

```
depths.zarr/
├── zarr.json           # Zarr v3 metadata
├── chr1/               # Per-chromosome arrays
│   └── ...
├── chr2/
│   └── ...
└── ...
```

**Root Metadata (`clam_metadata`):**

```json
{
  "contigs": [
    {"name": "chr1", "length": 248956422},
    {"name": "chr2", "length": 242193529}
  ],
  "column_names": ["sample1", "sample2", "sample3"],
  "chunk_size": 1000000
}
```

| Field | Type | Description |
|-------|------|-------------|
| `contigs` | array | List of chromosomes with names and lengths |
| `column_names` | array | Sample names in column order |
| `chunk_size` | integer | Chunk size in base pairs |

**Array Properties:**

| Property | Value |
|----------|-------|
| Data type | `uint32` |
| Shape | `(chromosome_length, num_samples)` |
| Chunk shape | `(chunk_size, num_samples)` |
| Compression | Blosc (Zstd, level 5, byte shuffle) |

**Per-Chromosome Metadata (`contig_info`):**

```json
{
  "contig": "chr1",
  "length": 248956422
}
```

---

## clam loci Output

### Callable Sites Zarr Store (default)

By default, `clam loci` produces a Zarr store containing callable sample counts per population.

**Structure:** Same as depth Zarr (see above).

**Root Metadata (`clam_metadata`):**

```json
{
  "contigs": [
    {"name": "chr1", "length": 248956422},
    {"name": "chr2", "length": 242193529}
  ],
  "column_names": ["PopA", "PopB", "PopC"],
  "chunk_size": 1000000
}
```

| Field | Type | Description |
|-------|------|-------------|
| `contigs` | array | List of chromosomes with names and lengths |
| `column_names` | array | Population names in column order |
| `chunk_size` | integer | Chunk size in base pairs |

**Array Properties:**

| Property | Value |
|----------|-------|
| Data type | `uint16` |
| Shape | `(chromosome_length, num_populations)` |
| Chunk shape | `(chunk_size, num_populations)` |
| Compression | Blosc (Zstd, level 5, byte shuffle) |

**Data interpretation:**

Each value represents the number of samples in that population that are callable at that position. For example, if `PopA` has 10 samples and the value at position 1000 is 8, then 8 of the 10 samples in `PopA` are callable at position 1000.

### Per-Sample Mask Zarr Store (`--per-sample`)

When using `--per-sample`, `clam loci` outputs a boolean mask indicating callability for each sample.

**Root Metadata (`clam_metadata`):**

```json
{
  "contigs": [
    {"name": "chr1", "length": 248956422}
  ],
  "column_names": ["sample1", "sample2", "sample3"],
  "chunk_size": 1000000
}
```

**Array Properties:**

| Property | Value |
|----------|-------|
| Data type | `bool` |
| Shape | `(chromosome_length, num_samples)` |
| Chunk shape | `(chunk_size, num_samples)` |
| Compression | PackBits |

**Data interpretation:**

- `true`: Site is callable for this sample
- `false`: Site is not callable for this sample

---

## clam stat Output

`clam stat` produces tab-separated (TSV) files in the specified output directory.

### pi.tsv

Nucleotide diversity (π) per population per window.

**Always generated.**

| Column | Type | Description |
|--------|------|-------------|
| `chrom` | string | Chromosome name |
| `start` | integer | Window start position (1-based) |
| `end` | integer | Window end position (1-based, inclusive) |
| `population` | string | Population name |
| `pi` | float | Nucleotide diversity estimate |
| `comparisons` | integer | Total pairwise comparisons (denominator) |
| `differences` | integer | Pairwise differences observed (numerator) |

**Example:**

```
chrom	start	end	population	pi	comparisons	differences
chr1	1	10000	PopA	0.00234	4500000	10530
chr1	1	10000	PopB	0.00198	3800000	7524
chr1	10001	20000	PopA	0.00241	4600000	11086
```

**Notes:**

- `pi = differences / comparisons`
- `comparisons` includes both variant and invariant callable sites
- `NaN` indicates no valid comparisons in the window

### pi_non_roh.tsv

Nucleotide diversity excluding samples in ROH regions (π~non-ROH~).

**Generated only when `--roh` is provided.**

Same columns as `pi.tsv`. At each site, samples that fall within an ROH region are excluded from the calculation. This provides a measure of non-ROH heterozygosity, which can serve as a proxy for the inbreeding load in a population ([Kardos et al. 2025](https://www.sciencedirect.com/science/article/pii/S016953472500182X)).

### dxy.tsv

Absolute divergence (d~xy~) between population pairs.

**Generated only when >1 population is defined.**

| Column | Type | Description |
|--------|------|-------------|
| `chrom` | string | Chromosome name |
| `start` | integer | Window start position (1-based) |
| `end` | integer | Window end position (1-based, inclusive) |
| `population1` | string | First population name |
| `population2` | string | Second population name |
| `dxy` | float | Absolute divergence estimate |
| `comparisons` | integer | Total between-population comparisons |
| `differences` | integer | Between-population differences |

**Example:**

```
chrom	start	end	population1	population2	dxy	comparisons	differences
chr1	1	10000	PopA	PopB	0.00312	8500000	26520
chr1	1	10000	PopA	PopC	0.00287	7200000	20664
chr1	1	10000	PopB	PopC	0.00301	6800000	20468
```

**Notes:**

- One row per population pair per window
- Population pairs are unordered (`PopA-PopB` is the same as `PopB-PopA`)
- `dxy = differences / comparisons`

### fst.tsv

Fixation index (F~ST~) between population pairs.

**Generated only when >1 population is defined.**

| Column | Type | Description |
|--------|------|-------------|
| `chrom` | string | Chromosome name |
| `start` | integer | Window start position (1-based) |
| `end` | integer | Window end position (1-based, inclusive) |
| `population1` | string | First population name |
| `population2` | string | Second population name |
| `fst` | float | F~ST~ estimate (Hudson estimator) |

**Example:**

```
chrom	start	end	population1	population2	fst
chr1	1	10000	PopA	PopB	0.156
chr1	1	10000	PopA	PopC	0.089
chr1	1	10000	PopB	PopC	0.201
```

**Notes:**

- F~ST~ is calculated using the Hudson estimator: F~ST~ = (d~xy~ - π~avg~) / d~xy~
- Values range from 0 (no differentiation) to 1 (complete differentiation)
- Negative values can occur when within-population diversity exceeds between-population divergence
- `NaN` indicates undefined F~ST~ (e.g., no polymorphic sites)

---

## Output Summary

| Command | Output | Format | Description |
|---------|--------|--------|-------------|
| `collect` | `*.zarr/` | Zarr | Raw depth values (uint32) |
| `loci` | `*.zarr/` | Zarr | Callable counts (uint16) or masks (bool) |
| `stat` | `pi.tsv` | TSV | Nucleotide diversity |
| `stat` | `pi_non_roh.tsv` | TSV | π excluding ROH (if `--roh`) |
| `stat` | `dxy.tsv` | TSV | Absolute divergence (if >1 pop) |
| `stat` | `fst.tsv` | TSV | Fixation index (if >1 pop) |

---

## Working with Outputs

### Reading Zarr in Python

```python
import zarr

# Open store
store = zarr.open("callable.zarr", mode="r")

# Get metadata
meta = store.attrs["clam_metadata"]
populations = meta["column_names"]
chunk_size = meta["chunk_size"]

# Read chromosome data
chr1 = store["chr1"][:]  # Shape: (chrom_length, num_populations)

# Read specific region (efficient)
region = store["chr1"][1000000:2000000, :]
```

For more detailed examples of working with Zarr outputs, including exploring depth distributions, see the [example notebooks](https://github.com/cademirch/clam/tree/main/notebooks) (coming soon).

### Reading TSV in Python

```python
import pandas as pd

pi = pd.read_csv("results/pi.tsv", sep="\t")
dxy = pd.read_csv("results/dxy.tsv", sep="\t")
fst = pd.read_csv("results/fst.tsv", sep="\t")
```

### Reading TSV in R

```r
pi <- read.table("results/pi.tsv", header=TRUE, sep="\t")
dxy <- read.table("results/dxy.tsv", header=TRUE, sep="\t")
fst <- read.table("results/fst.tsv", header=TRUE, sep="\t")
```
