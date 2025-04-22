---
title: Overview
---

# Overview

clam is a command-line tool for population genomics that provides an efficient workflow for calculating accurate population genetic statistics. Instead of requiring an all-sites VCF (which can be prohibitively large), clam uses a two-step process:

1. First, identify genomic regions with sufficient sequencing depth to be considered "callable" using `clam loci`
2. Then, calculate population genetic statistics using these callable regions and a VCF containing variants with `clam stat`

This approach eliminates the need for generating all-sites VCF files while still producing accurate diversity estimates, making it particularly suitable for large population genomics datasets.

