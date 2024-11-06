# Command-Line Help for `clam`

This document contains the help content for the `clam` command-line program.

**Command Overview:**

* [`clam`↴](#clam)
* [`clam loci`↴](#clam-loci)
* [`clam stat`↴](#clam-stat)

## `clam`

Callable Loci and More

**Usage:** `clam <COMMAND>`

###### **Subcommands:**

* `loci` — Calculate callable sites from depth statistics.
* `stat` — 



## `clam loci`

Calculate callable sites from depth statistics.

**Usage:** `clam loci [OPTIONS] <INFILE> <OUTPREFIX>`

###### **Arguments:**

* `<INFILE>` — Path to input D4 file
* `<OUTPREFIX>` — Output file prefix. The extension will be added automatically based on the `--no-counts` flag

###### **Options:**

* `-m`, `--min-depth <MIN_DEPTH>` — Minimum depth to consider site callable per individual

  Default value: `0`
* `-M`, `--max-depth <MAX_DEPTH>` — Maximum depth to consider site callable per individual

  Default value: `inf`
* `-d`, `--depth-proportion <DEPTH_PROPORTION>` — Proportion of samples passing thresholds at site to consider callable. Ignored when outputting counts

  Default value: `0`
* `-u`, `--min-mean-depth <MEAN_DEPTH_MIN>` — Minimum mean depth across all samples at site to consider callable. Ignored when outputting counts

  Default value: `0`
* `-U`, `--max-mean-depth <MEAN_DEPTH_MAX>` — Maximum mean depth across all samples at site to consider callable. Ignored when outputting counts

  Default value: `inf`
* `--no-counts` — Disable outputting counts; produces a .bed file instead

  Default value: `false`
* `-t`, `--threads <THREADS>` — Number of threads to use

  Default value: `1`
* `-p`, `--population-file <POPULATION_FILE>` — Path to file that defines populations. Tab separated: sample, population_name
* `--thresholds-file <THRESHOLD_FILE>` — Path to file that defines per-chromosome individual level thresholds. Tab separated: chrom, min, max
* `-x <EXCLUDE>` — Comma separated list of chromosomes to exclude
* `--exclude-file <EXCLUDE_FILE>` — Path to file with chromosomes to exclude, one per line



## `clam stat`

**Usage:** `clam stat <FILE>`

###### **Arguments:**

* `<FILE>`



<hr/>

<small><i>
    This document was generated automatically by
    <a href="https://crates.io/crates/clap-markdown"><code>clap-markdown</code></a>.
</i></small>

