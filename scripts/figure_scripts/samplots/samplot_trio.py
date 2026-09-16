#!/usr/bin/env python3
#
# Generate samplot trio images (proband/mother/father) for each unique SV in
# an R135 green-gene candidate TSV (as produced by filter_annotsv_panel.R).
#
# AnnotSV's "split" mode writes one row per gene a variant overlaps, so the
# same SV can appear on multiple rows here (e.g. an SV spanning two panel
# genes) - those are collapsed into a single plot, with every hit gene
# folded into the output filename/title.
#
# Usage:
#   samplot_trio.py \
#       --tsv <family_denovo_or_inherited_R135_green.tsv> \
#       --label denovo|inherited \
#       --family <FAMILY_ID> \
#       --proband-cram <path> --proband-id <SM-xxx> \
#       --mother-cram  <path> --mother-id  <SM-xxx> \
#       --father-cram  <path> --father-id  <SM-xxx> \
#       --reference <ref.fasta> \
#       --outdir <dir>

import argparse
import csv
import os
import subprocess
import sys


def safe(s):
    return "".join(c if c.isalnum() or c in "-_." else "_" for c in s)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--tsv", required=True)
    p.add_argument("--label", required=True)
    p.add_argument("--family", required=True)
    p.add_argument("--proband-cram", required=True)
    p.add_argument("--mother-cram", required=True)
    p.add_argument("--father-cram", required=True)
    p.add_argument("--proband-id", required=True)
    p.add_argument("--mother-id", required=True)
    p.add_argument("--father-id", required=True)
    p.add_argument("--reference", required=True)
    p.add_argument("--outdir", required=True)
    p.add_argument("--samplot-bin", default="samplot")
    args = p.parse_args()

    if not os.path.isfile(args.tsv):
        print(f"WARNING: {args.tsv} not found, skipping", file=sys.stderr)
        return

    with open(args.tsv) as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))

    if not rows:
        print(f"{args.tsv}: no variants, nothing to plot")
        return

    seen = {}
    order = []
    for row in rows:
        key = (row["SV_chrom"], row["SV_start"], row["SV_end"], row["SV_type"])
        if key not in seen:
            seen[key] = []
            order.append(key)
        gene = row["R135_green_genes"]
        if gene not in seen[key]:
            seen[key].append(gene)

    print(f"{args.tsv}: {len(rows)} rows -> {len(order)} unique SVs")
    os.makedirs(args.outdir, exist_ok=True)

    failures = 0
    for chrom, start, end, svtype in order:
        genes = "-".join(seen[(chrom, start, end, svtype)])
        chrom_ucsc = chrom if chrom.startswith("chr") else f"chr{chrom}"
        out_png = os.path.join(
            args.outdir,
            f"{args.family}_{args.label}_{safe(genes)}_{chrom_ucsc}_{start}_{end}_{svtype}.png",
        )
        cmd = [
            args.samplot_bin, "plot",
            "-n", f"Proband_{args.proband_id}", f"Mother_{args.mother_id}", f"Father_{args.father_id}",
            "-b", args.proband_cram, args.mother_cram, args.father_cram,
            "-r", args.reference,
            "-c", chrom_ucsc, "-s", start, "-e", end, "-t", svtype,
            "-o", out_png,
        ]
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            failures += 1
            last_line = result.stderr.strip().splitlines()[-1] if result.stderr.strip() else "unknown error"
            print(f"  FAILED {genes} {chrom_ucsc}:{start}-{end} ({svtype}): {last_line}", file=sys.stderr)
        else:
            print(f"  {genes} {chrom_ucsc}:{start}-{end} ({svtype}) -> {os.path.basename(out_png)}")

    print(f"{args.tsv}: {len(order) - failures}/{len(order)} plots written to {args.outdir}")


if __name__ == "__main__":
    main()
