---
title: Home
---

# clam

clam identifies genomic regions with sufficient sequencing depth to be considered "callable" and uses this information to calculate population genetic statistics from VCFs. It eliminates the need to generate all-sites VCF files while still producing accurate diversity estimates.

clam was designed specifically for large population genomics datasets.

## Installation

From bioconda:

```bash
conda create -n clam bioconda::clam
```

From source:

```bash
git clone https://github.com/cademirch/clam.git
cd clam
cargo build --release
./target/release/clam --help
```

## Quick Start

### 1. Generate callable loci from depth data

```bash
clam loci -o callable.zarr -m 10 -t 8 *.d4.gz
```

### 2. Calculate population genetic statistics

```bash
clam stat -o results/ -w 10000 -c callable.zarr variants.vcf.gz
```

## Commands

clam has three main commands:

| Command | Description |
|---------|-------------|
| `loci` | Generate callable loci from depth files (D4, GVCF) |
| `stat` | Calculate π, d~xy~, F~ST~ from VCF using callable sites |
| `collect` | Pre-process depth files into Zarr for repeated use |

## Citation

If you use clam in your research, please cite:

> Mirchandani C, Enbody E, Sackton TB, Corbett-Detig R. Efficient estimation of nucleotide diversity and divergence using callable loci (and more). *Mol Biol Evol*. 2025;msaf282. doi:[10.1093/molbev/msaf282](https://doi.org/10.1093/molbev/msaf282)

```bibtex
@article{Mirchandani2025-hx,
  title     = {Efficient estimation of nucleotide diversity and divergence using
               callable loci (and more)},
  author    = {Mirchandani, Cade and Enbody, Erik and Sackton, Timothy B and
               Corbett-Detig, Russ},
  journal   = {Mol. Biol. Evol.},
  publisher = {Oxford University Press (OUP)},
  number    = {msaf282},
  month     = nov,
  year      = 2025,
  language  = {en}
}
```

## License

clam is distributed under the terms of the [MIT](https://spdx.org/licenses/MIT.html) license.
