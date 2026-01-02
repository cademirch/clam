"""
Plot benchmarking results comparing clam collect vs d4tools merge.

Generates two plots:
1. runtime.pdf - Runtime vs sample count
2. size.pdf - Output size vs sample count
"""

import os
from pathlib import Path

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


def get_dir_size(path: Path) -> int:
    """Get total size of a directory in bytes."""
    total = 0
    for entry in path.rglob("*"):
        if entry.is_file():
            total += entry.stat().st_size
    return total


def get_file_size(path: Path) -> int:
    """Get file size in bytes."""
    return path.stat().st_size


def parse_benchmark_file(path: Path) -> pd.DataFrame:
    """Parse a Snakemake benchmark TSV file."""
    return pd.read_csv(path, sep="\t")


def collect_benchmark_data(
    sample_counts: list[int],
    clam_threads: list[int],
    benchmark_dir: Path,
) -> pd.DataFrame:
    """Collect runtime data from all benchmark files."""
    records = []

    # d4tools merge benchmarks
    for n in sample_counts:
        bench_file = benchmark_dir / "d4tools_merge" / f"n{n}.tsv"
        if bench_file.exists():
            df = parse_benchmark_file(bench_file)
            for _, row in df.iterrows():
                records.append(
                    {
                        "tool": "d4tools merge",
                        "samples": n,
                        "threads": 1,  # d4tools doesn't support threads
                        "time_s": row["s"],
                    }
                )

    # clam collect benchmarks
    for n in sample_counts:
        for t in clam_threads:
            bench_file = benchmark_dir / "clam_collect" / f"n{n}_t{t}.tsv"
            if bench_file.exists():
                df = parse_benchmark_file(bench_file)
                for _, row in df.iterrows():
                    records.append(
                        {
                            "tool": f"clam collect (t={t})",
                            "samples": n,
                            "threads": t,
                            "time_s": row["s"],
                        }
                    )

    return pd.DataFrame(records)


def collect_size_data(
    sample_counts: list[int],
    clam_threads: list[int],
    merged_d4_dir: Path,
    zarr_dir: Path,
) -> pd.DataFrame:
    """Collect output size data."""
    records = []

    # d4tools merge output sizes
    for n in sample_counts:
        d4_file = merged_d4_dir / f"n{n}.d4"
        if d4_file.exists():
            size_bytes = get_file_size(d4_file)
            records.append(
                {
                    "tool": "d4tools merge",
                    "samples": n,
                    "size_mb": size_bytes / (1024 * 1024),
                }
            )

    # clam collect output sizes (use t=1, size should be same for all thread counts)
    for n in sample_counts:
        zarr_path = zarr_dir / f"n{n}_t{clam_threads[0]}"
        if zarr_path.exists():
            size_bytes = get_dir_size(zarr_path)
            records.append(
                {
                    "tool": "clam collect",
                    "samples": n,
                    "size_mb": size_bytes / (1024 * 1024),
                }
            )

    return pd.DataFrame(records)


def plot_runtime(df: pd.DataFrame, output_path: Path) -> None:
    """Create runtime comparison plot."""
    sns.set_theme(style="whitegrid")
    fig, ax = plt.subplots(figsize=(10, 6))

    sns.lineplot(
        data=df,
        x="samples",
        y="time_s",
        hue="tool",
        style="tool",
        markers=True,
        dashes=False,
        errorbar="sd",
        ax=ax,
    )

    ax.set_xlabel("Number of Samples")
    ax.set_ylabel("Time (seconds)")
    ax.set_title("Runtime: clam collect vs d4tools merge")
    ax.legend(title="Tool", bbox_to_anchor=(1.02, 1), loc="upper left")

    plt.tight_layout()
    fig.savefig(output_path, bbox_inches="tight", dpi=150)
    plt.close(fig)


def plot_size(df: pd.DataFrame, output_path: Path) -> None:
    """Create output size comparison plot."""
    sns.set_theme(style="whitegrid")
    fig, ax = plt.subplots(figsize=(10, 6))

    sns.lineplot(
        data=df,
        x="samples",
        y="size_mb",
        hue="tool",
        style="tool",
        markers=True,
        dashes=False,
        ax=ax,
    )

    ax.set_xlabel("Number of Samples")
    ax.set_ylabel("Size (MB)")
    ax.set_title("Output Size: clam collect vs d4tools merge")
    ax.legend(title="Tool", bbox_to_anchor=(1.02, 1), loc="upper left")

    plt.tight_layout()
    fig.savefig(output_path, bbox_inches="tight", dpi=150)
    plt.close(fig)


def main() -> None:
    # Get parameters from Snakemake
    sample_counts = snakemake.params.sample_counts
    clam_threads = snakemake.params.clam_threads

    # Directories (relative to working directory)
    benchmark_dir = Path("benchmarks")
    merged_d4_dir = Path("results/merged_d4")
    zarr_dir = Path("results/zarr")

    # Output paths
    runtime_plot = Path(snakemake.output.runtime_plot)
    size_plot = Path(snakemake.output.size_plot)

    # Ensure output directory exists
    runtime_plot.parent.mkdir(parents=True, exist_ok=True)

    # Collect data
    runtime_df = collect_benchmark_data(sample_counts, clam_threads, benchmark_dir)
    size_df = collect_size_data(sample_counts, clam_threads, merged_d4_dir, zarr_dir)

    # Generate plots
    plot_runtime(runtime_df, runtime_plot)
    plot_size(size_df, size_plot)


if __name__ == "__main__":
    main()
else:
    # Called from Snakemake script directive
    main()
