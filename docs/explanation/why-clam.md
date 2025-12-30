---
title: Why clam?
---

# Why clam?

clam solves a fundamental problem in population genomics: how to accurately estimate genetic diversity when you don't have genotype calls at every site in the genome.

## The Problem: All-Sites VCFs Don't Scale

The gold standard for calculating population genetic statistics like nucleotide diversity (π), absolute divergence (d~xy~), and F~ST~ requires knowing the genotype at every site in the genome, not just variable sites. Traditionally, this meant generating an "all-sites" VCF that includes both variant and invariant positions.

For small datasets, this works fine. But for large population genomics projects with hundreds or thousands of samples, all-sites VCFs become impractical:

- **Storage**: An all-sites VCF for 100 samples across a 3Gb genome can easily exceed 1TB
- **Processing time**: Generating and parsing these files takes days or weeks
- **Memory**: Many tools can't handle files this large

## The Problem with Variants-Only VCFs

The obvious alternative is to use a variants-only VCF, which only contains positions where at least one sample differs from the reference. These files are much smaller and faster to work with.

But variants-only VCFs have a critical limitation: you can't distinguish between "this site is reference in all samples" and "this site has no data."

Consider a position where:

- Sample A has 30x coverage and is homozygous reference
- Sample B has 0x coverage (no data)

In a variants-only VCF, neither sample appears at this position. But when calculating diversity, these two situations are completely different:

- Sample A should contribute to the denominator (it's a valid comparison)
- Sample B should not (we don't know what allele it has)

Ignoring this distinction leads to systematic underestimation of diversity statistics, particularly in datasets with variable coverage or missing data ([Korunes and Samuk 2021](https://onlinelibrary.wiley.com/doi/full/10.1111/1755-0998.13326)).

## The Solution: Decouple Depth from Variants

clam takes a different approach: instead of storing genotypes at every site, it stores **callable site counts** separately from variants.

The workflow has two steps:

1. **`clam loci`**: Process depth information (from D4 files or GVCFs) to identify which sites are "callable" for each sample or population. A site is callable if it has sufficient sequencing depth to reliably call a genotype.

2. **`clam stat`**: Combine the callable site information with a variants-only VCF to calculate accurate population genetic statistics.

This approach gives you the accuracy of all-sites analysis with the efficiency of variants-only VCFs.

## Why This Works

The key insight is that for calculating diversity statistics, you don't need the actual genotype at invariant sites. You only need to know:

1. **At variant sites**: What are the allele counts?
2. **At invariant sites**: How many valid comparisons can be made?

The number of valid comparisons at invariant sites depends only on how many samples are callable, which clam tracks efficiently using a compact Zarr format.

For example, to calculate π (nucleotide diversity):

$$\pi = \frac{\text{number of pairwise differences}}{\text{number of pairwise comparisons}}$$

- The numerator comes from variant sites in the VCF
- The denominator comes from both variant sites AND callable invariant sites

By tracking callable sites separately, clam can compute the correct denominator without storing genotypes at millions of invariant positions.

## When to Use clam

clam is particularly useful when:

- You have many samples (tens to thousands)
- Your genome is large
- You have variable coverage across samples
- You want to run multiple analyses with different filtering criteria
- Storage or compute resources are limited

If you have a small dataset where all-sites VCFs are practical, you may not need clam. But for large-scale population genomics, clam provides a scalable path to accurate diversity estimates.
