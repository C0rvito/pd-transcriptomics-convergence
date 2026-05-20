#!/usr/bin/env python3
"""
merge_star_counts.py
--------------------
Merge STAR ReadsPerGene.out.tab files into one count matrix per gene.

Uses column 4 (stranded reverse strand), skipping the first 4 summary rows
produced by STAR (N_unmapped, N_multimapping, N_noFeature, N_ambiguous).

Expected input layout:
    <star_dir>/<GENE>/<SAMPLE>/<SAMPLE>.ReadsPerGene.out.tab

Output:
    <outdir>/<GENE>/star_counts_brutas.tsv
    (rows = Ensembl gene IDs, columns = sample IDs)

Usage:
    python bin/merge_star_counts.py
    python bin/merge_star_counts.py --star-dir data/star --outdir data/de/counts
"""

import argparse
import sys
from pathlib import Path

import pandas as pd


STAR_SUMMARY_ROWS = 4   # N_unmapped / N_multimapping / N_noFeature / N_ambiguous
COUNT_COL = 3           # 0-based index of column 4 (stranded reverse)


def collect_counts(gene_dir: Path) -> pd.DataFrame:
    """Read all ReadsPerGene.out.tab files for one gene and return a merged DataFrame."""
    frames = {}

    for sample_dir in sorted(gene_dir.iterdir()):
        if not sample_dir.is_dir():
            continue

        tab_files = list(sample_dir.glob("*.ReadsPerGene.out.tab"))
        if not tab_files:
            print(f"  [WARN] No ReadsPerGene.out.tab found in {sample_dir}, skipping.",
                  file=sys.stderr)
            continue

        tab_file = tab_files[0]
        sample_id = sample_dir.name

        df = pd.read_csv(
            tab_file,
            sep="\t",
            header=None,
            skiprows=STAR_SUMMARY_ROWS,
            usecols=[0, COUNT_COL],
            names=["gene_id", sample_id],
        )
        df = df.set_index("gene_id")
        frames[sample_id] = df[sample_id]

    if not frames:
        raise RuntimeError(f"No valid samples found in {gene_dir}")

    merged = pd.DataFrame(frames)
    merged.index.name = "gene_id"
    return merged


def main():
    parser = argparse.ArgumentParser(
        description="Merge STAR ReadsPerGene counts into one matrix per gene."
    )
    parser.add_argument(
        "--star-dir",
        type=Path,
        default=Path("data/star"),
        help="Root directory containing <GENE>/<SAMPLE>/ sub-folders (default: data/star)",
    )
    parser.add_argument(
        "--outdir",
        type=Path,
        default=Path("data/de/counts"),
        help="Output root directory; one sub-folder is created per gene (default: data/de/counts)",
    )
    args = parser.parse_args()

    if not args.star_dir.exists():
        sys.exit(f"[ERROR] star-dir not found: {args.star_dir}")

    gene_dirs = sorted(p for p in args.star_dir.iterdir() if p.is_dir())
    if not gene_dirs:
        sys.exit(f"[ERROR] No gene directories found under {args.star_dir}")

    for gene_dir in gene_dirs:
        gene = gene_dir.name
        print(f"Processing {gene}...")

        try:
            counts = collect_counts(gene_dir)
        except RuntimeError as exc:
            print(f"  [ERROR] {exc}", file=sys.stderr)
            continue

        out_dir = args.outdir / gene
        out_dir.mkdir(parents=True, exist_ok=True)
        out_file = out_dir / "star_counts_brutas.tsv"

        counts.to_csv(out_file, sep="\t")
        print(f"  -> {out_file}  ({counts.shape[0]} genes x {counts.shape[1]} samples)")

    print("Done.")


if __name__ == "__main__":
    main()
