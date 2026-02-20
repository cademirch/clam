---
title: Collect Depth Data
---

# Collect Depth Data

This guide covers how to use `clam collect` to pre-process depth files into a Zarr store for efficient reuse.

## When to Use collect

Use `clam collect` when you want to:

- **Run `loci` multiple times** with different threshold parameters
- **Explore depth distributions** across samples to decide on appropriate cutoffs
- **Save storage space** (Zarr is often smaller than merged D4 files)
- **Improve performance** for repeated analyses

If you only need to run `loci` once with known thresholds, you can skip `collect` and pass depth files directly to `loci`.

## Prerequisites

- Depth files in one of the supported formats:
    - Per-sample D4 files (uncompressed or bgzipped with index)
    - Merged D4 files (uncompressed only)
    - GVCF files (bgzipped and tabix indexed)

## Basic Usage

```bash
clam collect -o depths.zarr *.d4.gz
```

This reads all D4 files and stores the raw depth values in a Zarr store.

## Input Formats

### D4 Files

```bash
clam collect -o depths.zarr sample1.d4.gz sample2.d4.gz sample3.d4.gz
```

See [Input Formats](../reference/input-formats.md#d4-files) for details on D4 file requirements.

### Merged D4 Files

```bash
clam collect -o depths.zarr merged_samples.d4
```

clam automatically detects merged D4 files and extracts all sample names. Note that merged D4 files must be uncompressed.

### GVCF Files

```bash
clam collect -o depths.zarr --min-gq 20 *.g.vcf.gz
```

For GVCFs, you can filter by genotype quality (GQ) to exclude low-confidence depth values.

## Filtering Chromosomes

Exclude chromosomes you don't need to reduce storage:

```bash
# Exclude mitochondria and unplaced contigs
clam collect -o depths.zarr -x chrM,chrUn *.d4.gz

# Only include autosomes
clam collect -o depths.zarr -i chr1,chr2,chr3,chr4,chr5 *.d4.gz
```

## Performance

Use multiple threads for faster processing:

```bash
clam collect -o depths.zarr -t 16 *.d4.gz
```

Adjust chunk size if needed (default 1Mb):

```bash
clam collect -o depths.zarr --chunk-size 500000 *.d4.gz
```

## Workflow Example

### Step 1: Collect Depth Once

```bash
clam collect -o depths.zarr -t 16 *.d4.gz
```

### Step 2: Explore Depth Distributions

Use Python to explore the depth distribution and decide on thresholds. See the [example notebooks](https://github.com/cademirch/clam/tree/main/notebooks) for detailed examples (coming soon).

```python
import zarr
import numpy as np

store = zarr.open("depths.zarr", mode="r")
chr1_depths = store["chr1"][:]

# Get per-sample mean depths
mean_depths = np.mean(chr1_depths, axis=0)
print(f"Mean depths per sample: {mean_depths}")
```

### Step 3: Run loci with Different Thresholds

```bash
# Conservative thresholds
clam loci -o callable_strict.zarr -m 15 -M 80 depths.zarr

# Relaxed thresholds
clam loci -o callable_relaxed.zarr -m 5 -M 150 depths.zarr

# Test different proportion requirements
clam loci -o callable_d50.zarr -m 10 -d 0.5 depths.zarr
clam loci -o callable_d80.zarr -m 10 -d 0.8 depths.zarr
```

!!! tip "Population Flow"
    If you specified populations during `collect` (via `--samples`), they are stored in the Zarr metadata. `clam loci` reads them automatically, so you don't need to re-specify `-p` or `--samples`:
    
    ```bash
    # Populations were defined during collect
    clam collect -o depths.zarr --samples samples.tsv
    
    # loci reads populations from zarr metadata — no -p needed
    clam loci -o callable.zarr -m 10 depths.zarr
    
    # Override with different populations if needed
    clam loci -o callable.zarr -m 10 depths.zarr --samples different_pops.tsv
    ```

### Step 4: Compare Results

```bash
# Run stat with each callable set
clam stat -o results_strict/ -w 10000 -c callable_strict.zarr variants.vcf.gz
clam stat -o results_relaxed/ -w 10000 -c callable_relaxed.zarr variants.vcf.gz
```

## Storage Efficiency

The Zarr format with Zstd compression typically achieves better compression than merged D4 files, especially for large sample counts. The exact savings depend on your data, but 20-50% reduction is common.

## Output Format

The output is a Zarr store containing:

- Raw depth values (uint32) for each sample at each position
- Metadata including sample names, chromosome lengths, chunk size, and population assignments (when specified)

See [Output Formats](../reference/output-formats.md#depth-zarr-store) for details on the Zarr structure.

## Complete Example

```bash
# Collect depth from 100 samples, excluding non-standard chromosomes
clam collect \
  -o depths.zarr \
  --exclude-file exclude_chroms.txt \
  --chunk-size 1000000 \
  -t 16 \
  *.d4.gz

# Now use the Zarr store for multiple loci runs
clam loci -o callable_m10.zarr -m 10 -M 100 -p pops.tsv depths.zarr
clam loci -o callable_m15.zarr -m 15 -M 100 -p pops.tsv depths.zarr
```

## Next Steps

After collecting depth data, use `clam loci` to generate callable sites:

- [Generate Callable Loci](generate-callable-loci.md)

Then calculate statistics with `clam stat`:

- [Calculate Statistics](calculate-statistics.md)
