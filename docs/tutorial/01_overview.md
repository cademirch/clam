---
title: Getting Started
---

# Getting Started

!!! note "Work in Progress"
    This tutorial is under development. Check back soon for a complete walkthrough using example data.

This tutorial will guide you through a complete clam workflow using example data, from raw depth files to population genetic statistics.

## What You'll Learn

- How to prepare input files for clam
- Generate callable loci with `clam loci`
- Calculate π, d~xy~, and F~ST~ with `clam stat`
- Interpret and visualize the results

## Prerequisites

- clam installed (see [Installation](../index.md#installation))
- Basic familiarity with command-line tools
- Python or R for visualization (optional)

## Example Dataset

This tutorial uses simulated data from the clam test suite, representing:

- 20 samples across 2 populations (Pop1, Pop2)
- 10 samples per population
- D4 depth files and a variants-only VCF

## Tutorial Outline

### Part 1: Prepare Your Data

- Organize depth files (D4 or GVCF)
- Create a population file
- Inspect chromosome names and lengths

### Part 2: Generate Callable Loci

- Choose appropriate depth thresholds
- Run `clam loci` with population definitions
- Examine the output Zarr store

### Part 3: Calculate Statistics

- Run `clam stat` with the callable loci
- Understand window-based calculations
- Explore the output TSV files

### Part 4: Analyze Results

- Load results in Python/R
- Plot diversity across the genome
- Compare populations

## Coming Soon

Full tutorial content with step-by-step instructions and example code.

In the meantime, see:

- [How-to: Generate Callable Loci](../how-to/generate-callable-loci.md)
- [How-to: Calculate Statistics](../how-to/calculate-statistics.md)
- [Core Concepts](../explanation/concepts.md)
