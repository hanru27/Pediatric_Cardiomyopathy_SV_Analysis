#!/usr/bin/env python3
"""
Annotates a candidate SV VCF with SpliceAI (breakpoint splice-impact) and
GREEN-DB (regulatory region overlap) scores, one TSV row per SV.

SpliceAI can't score Manta's symbolic ALT alleles (<DEL>, <INS>, ...)
directly, so each SV's two breakpoints are scored as minimal single-base
proxy deletions instead, keeping the max delta score across both. GREEN-DB
(via GREEN-VARAN's `sv` mode) runs against the real SV VCF directly.

Usage:
  annotate_splice_greendb.py --vcf candidates.vcf.gz --ref genome.fasta \
    --greenvaran-bin tools/GREEN-VARAN/bin/greenvaran \
    --greendb tools/GREEN-VARAN/resources/GRCh38/GRCh38_GREEN-DB.bed.gz \
    --greendb-schema tools/GREEN-VARAN/config/greendb_schema_v2.5.json \
    --spliceai-annotation grch38 \
    --out family_splice_greendb.tsv
"""
import argparse
import gzip
import os
import subprocess
import sys
import tempfile
from pathlib import Path

import pysam


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--vcf", required=True, help="candidate SV VCF (bgzipped, indexed)")
    p.add_argument("--ref", required=True, help="reference genome fasta")
    p.add_argument("--greenvaran-bin", required=True)
    p.add_argument("--greendb", required=True)
    p.add_argument("--greendb-schema", required=True)
    p.add_argument("--spliceai-annotation", default="grch38")
    p.add_argument("--spliceai-distance", type=int, default=50)
    p.add_argument("--out", required=True)
    return p.parse_args()


def sv_records(vcf_path):
    """Yield (sv_id, chrom, pos, end) for every record in the candidate VCF."""
    vf = pysam.VariantFile(vcf_path)
    for rec in vf:
        end = rec.info.get("END", rec.pos)
        yield rec.id, rec.chrom, rec.pos, int(end)


def build_breakpoint_proxy_vcf(vcf_path, ref_path, out_path):
    """One minimal 2bp->1bp deletion proxy variant per SV breakpoint (POS and END)."""
    fasta = pysam.FastaFile(ref_path)
    breakpoints = []  # (sv_id, which, chrom, pos)
    for sv_id, chrom, pos, end in sv_records(vcf_path):
        breakpoints.append((sv_id, "start", chrom, pos))
        if end != pos:
            breakpoints.append((sv_id, "end", chrom, end))

    with open(out_path, "w") as f:
        f.write("##fileformat=VCFv4.2\n")
        for name, length in zip(fasta.references, fasta.lengths):
            f.write(f"##contig=<ID={name},length={length}>\n")
        f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for sv_id, which, chrom, pos in sorted(breakpoints, key=lambda b: (b[2], b[3])):
            ref_seq = fasta.fetch(chrom, pos - 1, pos + 1).upper()  # 2 bases: pos, pos+1
            if len(ref_seq) < 2 or "N" in ref_seq:
                continue
            proxy_id = f"{sv_id}__{which}"
            f.write(f"{chrom}\t{pos}\t{proxy_id}\t{ref_seq}\t{ref_seq[0]}\t.\t.\t.\n")
    fasta.close()
    return out_path


def run_spliceai(proxy_vcf, ref_path, annotation, distance, out_vcf):
    # invoked via a full python interpreter path (not `conda activate`), so
    # PATH doesn't include this env's bin/ - resolve the spliceai console
    # script relative to the running interpreter instead of relying on PATH
    spliceai_bin = os.path.join(os.path.dirname(sys.executable), "spliceai")
    cmd = [
        spliceai_bin, "-I", proxy_vcf, "-O", out_vcf,
        "-R", ref_path, "-A", annotation, "-D", str(distance),
    ]
    subprocess.run(cmd, check=True)


def parse_spliceai_output(spliceai_vcf):
    """sv_id -> {"score": max delta score across both breakpoints, "gene": SYMBOL, "type": which delta was max}"""
    results = {}
    delta_names = ["DS_AG", "DS_AL", "DS_DG", "DS_DL"]
    with open(spliceai_vcf) as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            proxy_id, info = fields[2], fields[7]
            sv_id, which = proxy_id.rsplit("__", 1)
            sai = next((kv[10:] for kv in info.split(";") if kv.startswith("SpliceAI=")), None)
            if not sai or sai == ".":
                continue
            parts = sai.split("|")
            gene = parts[1]
            scores = [float(x) for x in parts[2:6]]
            max_idx = max(range(4), key=lambda i: scores[i])
            max_score = scores[max_idx]
            prev = results.get(sv_id)
            if prev is None or max_score > prev["score"]:
                results[sv_id] = {"score": max_score, "gene": gene, "type": delta_names[max_idx], "breakpoint": which}
    return results


def run_greenvaran(vcf_path, greenvaran_bin, greendb, greendb_schema, out_vcf):
    cmd = [
        greenvaran_bin, "sv",
        "-i", vcf_path, "-o", out_vcf,
        "-d", greendb, "-s", greendb_schema,
        "-u",
    ]
    subprocess.run(cmd, check=True)


def parse_greenvaran_output(greenvaran_vcf):
    results = {}
    opener = gzip.open if greenvaran_vcf.endswith(".gz") else open
    with opener(greenvaran_vcf, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            sv_id, info = fields[2], fields[7]
            info_kv = dict(kv.split("=", 1) for kv in info.split(";") if kv.startswith("greendb_") and "=" in kv)
            if not info_kv:
                continue
            results[sv_id] = {
                "id": info_kv.get("greendb_id", ""),
                "type": info_kv.get("greendb_stdtype", ""),
                "source": info_kv.get("greendb_dbsource", ""),
                "genes": info_kv.get("greendb_genes", ""),
                "constraint": info_kv.get("greendb_constraint", ""),
            }
    return results


def main():
    args = parse_args()

    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)

        proxy_vcf = build_breakpoint_proxy_vcf(args.vcf, args.ref, tmp / "proxy.vcf")
        spliceai_out = tmp / "proxy.spliceai.vcf"
        run_spliceai(str(proxy_vcf), args.ref, args.spliceai_annotation, args.spliceai_distance, str(spliceai_out))
        splice_results = parse_spliceai_output(str(spliceai_out))

        greenvaran_out = tmp / "greendb.vcf"
        run_greenvaran(args.vcf, args.greenvaran_bin, args.greendb, args.greendb_schema, str(greenvaran_out))
        greendb_results = parse_greenvaran_output(str(greenvaran_out))

    all_sv_ids = [sv_id for sv_id, *_ in sv_records(args.vcf)]

    with open(args.out, "w") as out:
        out.write("SV_ID\tSpliceAI_max_score\tSpliceAI_gene\tSpliceAI_effect\tSpliceAI_breakpoint\t"
                   "GREENDB_id\tGREENDB_type\tGREENDB_source\tGREENDB_genes\tGREENDB_constraint\n")
        for sv_id in all_sv_ids:
            sp = splice_results.get(sv_id, {})
            gd = greendb_results.get(sv_id, {})
            out.write("\t".join([
                sv_id,
                str(sp.get("score", "")),
                sp.get("gene", ""),
                sp.get("type", ""),
                sp.get("breakpoint", ""),
                gd.get("id", ""),
                gd.get("type", ""),
                gd.get("source", ""),
                gd.get("genes", ""),
                str(gd.get("constraint", "")),
            ]) + "\n")

    print(f"Wrote {len(all_sv_ids)} SVs -> {args.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
