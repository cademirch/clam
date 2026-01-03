---
title: Why clam?
---

# Why clam?

Accurate estimation of population genetic statistics like nucleotide diversity (π), divergence (d~xy~), and F~ST~ requires knowing which sites in the genome were actually assayed. For example, π is calculated as:

$$\pi = \frac{\text{number of pairwise differences}}{\text{number of pairwise comparisons}}$$

The denominator requires knowing how many sites were successfully genotyped. Getting this wrong leads to biased estimates.

## The Problem with Variants-Only VCFs

Standard VCF files only record positions where at least one sample differs from the reference. They don't record genotypes at invariant sites, which creates a problem: you can't distinguish between "this site is reference in all samples" and "this site has no data."

Consider a position where:

- Sample A has 30x coverage and is homozygous reference
- Sample B has 0x coverage (no data)

In a variants-only VCF, neither sample appears at this position. But when calculating diversity, these situations are different:

- Sample A should contribute to the denominator (it's a valid comparison)
- Sample B should not (we don't know what allele it has)

Many tools assume that sites absent from the VCF are homozygous reference, leading to systematic underestimation of diversity statistics ([Korunes and Samuk 2021](https://onlinelibrary.wiley.com/doi/full/10.1111/1755-0998.13326)). The estimators themselves aren't statistically biased; the apparent bias arises from mismeasurement of the denominator when incorrect assumptions about unobserved sites result in erroneous counts of total base pairs assayed.

## Why All-Sites VCFs Don't Scale

One solution is to generate an "all-sites" VCF that includes both variant and invariant positions. This correctly records which sites were successfully genotyped across samples, allowing proper enumeration of the denominator.

However, all-sites VCFs become impractical as datasets grow. For large population genomics projects with hundreds or thousands of samples, they pose computational challenges in creation, storage, and analysis. Additionally, filtering approaches for all-sites VCFs often differ from those applied to variants-only VCFs, introducing potentially unaccounted for confounding effects.

## clam's Approach: Callable Loci

clam uses callable loci—determined from sequencing depth—to estimate population genetic statistics from a variants-only VCF. The key insight is that for calculating diversity statistics, you don't need the actual genotype at invariant sites. You only need to know:

1. **At variant sites**: What are the allele counts?
2. **At invariant sites**: How many samples had sufficient coverage to call a genotype?

The workflow has two steps:

1. **`clam loci`**: Process depth information (from D4 files or GVCFs) to identify which sites are "callable" for each sample or population. A site is callable if it has sufficient sequencing depth to reliably call a genotype.

2. **`clam stat`**: Combine the callable site information with a variants-only VCF to calculate population genetic statistics.

By tracking callable sites separately, clam computes the correct denominator without storing genotypes at millions of invariant positions. The numerator comes from variant sites in the VCF; the denominator comes from both variant sites and callable invariant sites.

## When to Use clam

clam is useful when:

- You have many samples (tens to thousands)
- Your genome is large
- You have variable coverage across samples
- You want to run multiple analyses with different filtering criteria
- Storage or compute resources are limited

If you have a small dataset where all-sites VCFs are practical, you may not need clam. But for large-scale population genomics, clam provides a scalable path to accurate diversity estimates.
