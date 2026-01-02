# Benchmarking: clam collect vs d4tools merge

Compares `clam collect` and `d4tools merge` for combining per-sample depth files, measuring runtime and output size as sample count increases.

## Setup

Simulates 10 samples of paired-end reads (~10x coverage) from E. coli, converts to D4 format, then benchmarks merging 2-10 samples.

## Usage

```bash
cd benches
pixi run snakemake -s workflow/Snakefile --cores 8 --resources io_heavy=1
```

The `io_heavy=1` resource ensures benchmarked rules run sequentially to avoid I/O contention.

## Results

### Runtime

![Runtime comparison](results/plots/runtime.png)

### Output Size

![Size comparison](results/plots/size.png)
