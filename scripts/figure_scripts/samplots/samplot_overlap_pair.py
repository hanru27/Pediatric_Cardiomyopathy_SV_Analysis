#!/usr/bin/env python3
#
# Generate samplot images for the genome-wide overlap candidates (see
# find_overlapping_denovo.R / build_genewide_overlap_candidate_table.py):
# cross-family pairs of de novo SVs at the exact same locus, matching
# proband phenotype, flagged as possible recurrent artefacts/hotspots.
#
# Unlike samplot_trio.py (one family's proband+mother+father at one locus),
# this plots BOTH probands from a pair TOGETHER in a single samplot image -
# the useful comparison here is proband-vs-proband at the shared locus, to
# see by eye whether the two unrelated calls look like the same mapping/
# repeat artefact (near-identical pileup pattern) or genuinely independent
# events. Parents are omitted; de novo status for both is already confirmed
# upstream.
#
# Usage:
#   samplot_overlap_pair.py \
#       --candidates-tsv <de_novo_genewide_overlap_candidates.tsv> \
#       --trios <trios_all37.txt> \
#       --reference <ref.fasta> \
#       --outdir <dir> \
#       --samplot-bin <path/to/samplot>

import argparse
import csv
import os
import subprocess
import sys


def safe(s):
    return "".join(c if c.isalnum() or c in "-_." else "_" for c in s)


def load_trios(path):
    """Family -> dict(proband_id, proband_sm, proband_cram, ...)."""
    trios = {}
    with open(path, encoding="utf-8") as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 10:
                continue
            fam = parts[0]
            trios[fam] = {
                "proband_id": parts[1], "proband_sm": parts[2], "proband_cram": parts[3],
                "mother_id": parts[4], "mother_sm": parts[5], "mother_cram": parts[6],
                "father_id": parts[7], "father_sm": parts[8], "father_cram": parts[9],
            }
    return trios


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--candidates-tsv", required=True,
                    help="TSV with columns Gene_name,Family_A,Family_B,Chrom,Start,End,Type")
    p.add_argument("--trios", required=True)
    p.add_argument("--reference", required=True)
    p.add_argument("--outdir", required=True)
    p.add_argument("--samplot-bin", default="samplot")
    args = p.parse_args()

    trios = load_trios(args.trios)
    os.makedirs(args.outdir, exist_ok=True)

    with open(args.candidates_tsv, encoding="utf-8") as f:
        rows = list(csv.DictReader(f, delimiter="\t"))

    print(f"{args.candidates_tsv}: {len(rows)} candidate pairs")

    failures = 0
    for row in rows:
        gene = row["Gene_name"]
        fam_a, fam_b = row["Family_A"], row["Family_B"]
        chrom, start, end, svtype = row["Chrom"], row["Start"], row["End"], row["Type"]

        if fam_a not in trios or fam_b not in trios:
            missing = fam_a if fam_a not in trios else fam_b
            print(f"  SKIPPED {gene} {fam_a}/{fam_b}: family {missing} not in {args.trios}", file=sys.stderr)
            failures += 1
            continue

        a, b = trios[fam_a], trios[fam_b]
        chrom_ucsc = chrom if chrom.startswith("chr") else f"chr{chrom}"
        out_png = os.path.join(
            args.outdir,
            f"{safe(gene)}_{fam_a}_vs_{fam_b}_{chrom_ucsc}_{start}_{end}_{svtype}.png",
        )
        cmd = [
            args.samplot_bin, "plot",
            "-n", f"Proband_{fam_a}_{a['proband_sm']}", f"Proband_{fam_b}_{b['proband_sm']}",
            "-b", a["proband_cram"], b["proband_cram"],
            "-r", args.reference,
            "-c", chrom_ucsc, "-s", start, "-e", end, "-t", svtype,
            "-o", out_png,
        ]
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            failures += 1
            last_line = result.stderr.strip().splitlines()[-1] if result.stderr.strip() else "unknown error"
            print(f"  FAILED {gene} {fam_a}_vs_{fam_b} {chrom_ucsc}:{start}-{end} ({svtype}): {last_line}", file=sys.stderr)
        else:
            print(f"  {gene} {fam_a}_vs_{fam_b} {chrom_ucsc}:{start}-{end} ({svtype}) -> {os.path.basename(out_png)}")

    print(f"{len(rows) - failures}/{len(rows)} plots written to {args.outdir}")


if __name__ == "__main__":
    main()
