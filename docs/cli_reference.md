---
title: CLI Reference
---
# Command-Line Help for `clam`

This page contains the help content for the `clam` command-line program.

**Command Overview:**

* [`clam`↴](#clam)
* [`clam loci`↴](#clam-loci)
* [`clam stat`↴](#clam-stat)

## `clam`

```console
Usage: clam [OPTIONS] <COMMAND>

Commands:
  loci  Calculate callable sites from depth statistics.
  stat  Calculate population genetic statistics from VCF using callable sites.
  help  Print this message or the help of the given subcommand(s)

Options:
  -v, --verbose...  Increase verbosity (-v, -vv for more verbosity)
  -q, --quiet       Suppress output (overrides verbosity)
  -h, --help        Print help
  -V, --version     Print version
```

## `clam loci`

```console
Usage: clam loci [OPTIONS] -o <OUTDIR> [INPUT]...

Arguments:
  [INPUT]...  Input files (D4 format by default). Specify one or more files directly

Options:
  -f, --filelist <FILELIST>  Path to file containing list of input files, one per line. Use this instead of positional input arguments for many files
      --gvcf                 Use GVCF format instead of default D4 format for input files
      --merged               Input is a merged D4 file (single file containing multiple samples)
  -o <OUTDIR>                Output directory for results (required)
      --bed                  Write additional BED file. Note: This can be slow for large datasets
  -h, --help                 Print help
  -V, --version              Print version

Sample-level Thresholds:
  -m, --min-depth <MIN_DEPTH>
          Minimum depth to consider a site callable for each individual [default: 0]
  -M, --max-depth <MAX_DEPTH>
          Maximum depth to consider a site callable for each individual [default: inf]
      --thresholds-file <THRESHOLD_FILE>
          Custom thresholds per chromosome. Tab-separated file: chrom, min, max

Population-level Thresholds:
  -d, --depth-proportion <DEPTH_PROPORTION>
          Proportion of samples that must pass thresholds at a site to consider it callable. Value between 0.0 and 1.0 [default: 0]
  -u, --min-mean-depth <MEAN_DEPTH_MIN>
          Minimum mean depth across all samples required at a site to consider it callable [default: 0]
  -U, --max-mean-depth <MEAN_DEPTH_MAX>
          Maximum mean depth across all samples allowed at a site to consider it callable [default: inf]
  -p, --population-file <POPULATION_FILE>
          Path to file that defines populations. Tab separated: sample, population_name

Chromosome Filtering:
  -x <EXCLUDE>...                    Comma separated list of chromosomes to exclude. Example: --exclude chr1,chr2,chrX
      --exclude-file <EXCLUDE_FILE>  Path to file with chromosomes to exclude, one per line
  -i <INCLUDE>...                    Comma separated list of chromosomes to include (restrict analysis to). Example: --include chr1,chr2,chr3
      --include-file <INCLUDE_FILE>  Path to file with chromosomes to include, one per line

Performance:
  -t, --threads <THREADS>  Number of threads to use for parallel processing [default: 1]

EXAMPLES:
    # Basic usage with positional input files
    clam loci -o output_dir input1.d4 input2.d4
    
    # Using a file list instead of positional arguments
    clam loci -f filelist.txt -o output_dir
    
    # Set custom depth thresholds
    clam loci -o output_dir -m 10 -M 100 input1.d4 input2.d4
```
## `clam stat`

```console
Usage: clam stat [OPTIONS] <--window-size <WINDOW_SIZE>|--regions-file <REGIONS_FILE>> <VCF> [CALLABLE_SITES]

Arguments:
  <VCF>             Path to input VCF file
  [CALLABLE_SITES]  Path to input callable sites D4 file from clam loci

Options:
  -o, --outdir <OUTDIR>
          Where to write output files. Defaults to current working directory
  -t, --threads <THREADS>
          Number of threads to use [default: 1]
  -w, --window-size <WINDOW_SIZE>
          Size of windows for statistics in bp. Conflicts with 'regions-file'
  -r, --regions-file <REGIONS_FILE>
          File specifying regions to calculate statistics for. Conflicts with 'window-size'
  -s, --sites-file <SITES_FILE>
          Specify sites to consider for calculations. Bed format
  -p, --population-file <POPULATION_FILE>
          Path to file that defines populations. Tab separated: sample, population_name
  -f, --fai <FASTA_INDEX>
          Path to fasta index for reference VCF was called against. Only needed if VCF does not have contig info in the header
  -x <EXCLUDE>...
          Comma separated list of chromosomes to exclude
      --exclude-file <EXCLUDE_FILE>
          Path to file with chromosomes to exclude, one per line
      --roh-file <ROH_FILE>
          Path to RoH file
  -i <INCLUDE>...
          Comma separated list of chromosomes to include (restrict analysis to)
      --include-file <INCLUDE_FILE>
          Path to file with chromosomes to include, one per line
  -h, --help
          Print help
```

