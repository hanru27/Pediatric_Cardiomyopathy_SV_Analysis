#!/usr/bin/env python3
"""
Filters Overlapping_SV_pairs to Possible_artefact == "Yes" pairs, then
enriches each with its AnnotSV/SpliceAI/GREEN-DB annotation by Manta SV ID.

Usage:
  python3 build_genewide_overlap_candidate_table.py \
    --overlap-xlsx <denovo_overlap_analysis_AF_filtered.xlsx> \
    --denovo-calls-dir <results/De_novo_calls> \
    --out <de_novo_genewide_overlap_candidates.xlsx>
"""
import argparse
import csv
import sys
import zipfile
import re
import html as htmllib
from pathlib import Path

from openpyxl import Workbook
from openpyxl.styles import Font, PatternFill
from openpyxl.utils import get_column_letter


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--overlap-xlsx", required=True)
    p.add_argument("--denovo-calls-dir", required=True)
    p.add_argument("--out", required=True)
    return p.parse_args()


#xlsx reader 
def read_sheet(xlsx_path, sheet_name):
    with zipfile.ZipFile(xlsx_path) as z:
        wb_xml = z.read("xl/workbook.xml").decode("utf-8")
        sheets = re.findall(r'<sheet name="([^"]*)"[^/]*r:id="(rId\d+)"', wb_xml)
        rels_xml = z.read("xl/_rels/workbook.xml.rels").decode("utf-8")
        rels = dict(re.findall(r'Id="(rId\d+)"[^>]*Target="([^"]*)"', rels_xml))
        target = next(rels[rid] for name, rid in sheets if name == sheet_name)
        sheet_path = "xl/" + target if not target.startswith("xl/") else target

        shared = []
        if "xl/sharedStrings.xml" in z.namelist():
            ss_xml = z.read("xl/sharedStrings.xml").decode("utf-8")
            for si in re.findall(r"<si>(.*?)</si>", ss_xml, re.S):
                texts = re.findall(r"<t[^>]*>(.*?)</t>", si, re.S)
                shared.append("".join(htmllib.unescape(t) for t in texts))

        data = z.read(sheet_path).decode("utf-8")

    def col_idx(col):
        idx = 0
        for c in col:
            idx = idx * 26 + (ord(c) - ord("A") + 1)
        return idx - 1

    rows_out = {}
    for rnum, content in re.findall(r'<row r="(\d+)"[^>]*>(.*?)</row>', data, re.S):
        rowd = {}
        for col, _rr, attrs, inner in re.findall(
            r'<c r="([A-Z]+)(\d+)"([^>]*?)(?:/>|>(.*?)</c>)', content, re.S
        ):
            t = "s" if 't="s"' in attrs else ("b" if 't="b"' in attrs else "n")
            val = ""
            if inner:
                m = re.search(r"<v>(.*?)</v>", inner, re.S)
                if m:
                    v = m.group(1)
                    val = shared[int(v)] if t == "s" else ("TRUE" if (t == "b" and v == "1") else v)
            rowd[col_idx(col)] = val
        rows_out[int(rnum)] = rowd

    maxcol = max((max(r.keys()) for r in rows_out.values() if r), default=-1) + 1
    ordered = [rows_out[r] for r in sorted(rows_out)]
    table = [[row.get(c, "") for c in range(maxcol)] for row in ordered]
    header, body = table[0], table[1:]
    return header, [r for r in body if any(c != "" for c in r)]


def load_family_tsv(path):
    with open(path, encoding="utf-8") as f:
        r = csv.reader(f, delimiter="\t")
        header = next(r)
        rows = list(r)
    return {h: i for i, h in enumerate(header)}, rows


def load_splice_greendb(path):
    idx, rows = load_family_tsv(path)
    return idx, {row[idx["SV_ID"]]: row for row in rows}


def find_annotsv_row(fam_idx, fam_rows, gene, chrom, start, end):
    for row in fam_rows:
        if (
            row[fam_idx["Gene_name"]] == gene
            and row[fam_idx["SV_chrom"]] == chrom
            and row[fam_idx["SV_start"]] == start
            and row[fam_idx["SV_end"]] == end
        ):
            return row
    return None


OUT_COLS = [
    "Gene_name", "Family_A", "Family_B", "Chrom", "Start", "End", "Type", "GT",
    "Reciprocal_pct_overlap", "SV_length", "Location", "Tandem",
    "HI", "TS", "DDD_HI_percent", "GnomAD_pLI", "LOEUF_bin",
    "Overlapped_CDS_percent", "Frameshift", "Dist_nearest_SS", "Nearest_SS_type",
    "SpliceAI_max_score", "GREENDB_type", "GREENDB_constraint", "OMIM_phenotype",
    "N_families_exact_locus", "Phenotype_A", "Phenotype_B", "Artefact_reason",
]


def main():
    args = parse_args()
    header, rows = read_sheet(args.overlap_xlsx, "Overlapping_SV_pairs")
    idx = {h: i for i, h in enumerate(header)}

    flagged = [r for r in rows if r[idx["Possible_artefact"]] == "Yes"]
    print(f"{len(flagged)}/{len(rows)} same-phenotype overlapping pairs flagged "
          f"(exact 100% breakpoints, >=3-family recurrence)", file=sys.stderr)

    calls_dir = Path(args.denovo_calls_dir)
    fam_cache, splice_cache = {}, {}

    def annotsv(fam):
        if fam not in fam_cache:
            fam_cache[fam] = load_family_tsv(
                calls_dir / fam / "results/annotsv/denovo_allgenes_AF_filtered.tsv"
            )
        return fam_cache[fam]

    def splice(fam):
        if fam not in splice_cache:
            splice_cache[fam] = load_splice_greendb(
                calls_dir / fam / "results/annotsv/denovo_splice_greendb.tsv"
            )
        return splice_cache[fam]

    out_rows = []
    for r in flagged:
        gene = r[idx["Gene_name"]]
        famA, famB = r[idx["Family_A"]], r[idx["Family_B"]]
        chrom, start, end = r[idx["Chrom_A"]], r[idx["Start_A"]], r[idx["End_A"]]
        typ, gt = r[idx["Type_A"]], r[idx["GT_A"]]
        recip = r[idx["Reciprocal_pct_overlap"]]
        phA, phB = r[idx["Phenotype_A"]], r[idx["Phenotype_B"]]
        reason = r[idx["Artefact_reason"]]
        nfam = reason.split("shared by ")[1].split(" unrelated")[0] if "shared by" in reason else ""

        fidx, frows = annotsv(famA)
        arow = find_annotsv_row(fidx, frows, gene, chrom, start, end)
        if arow is None:
            print(f"WARNING: no AnnotSV row for {gene} {famA} chr{chrom}:{start}-{end}", file=sys.stderr)
            continue

        sidx, sdict = splice(famA)
        srow = sdict.get(arow[fidx["ID"]])

        tandem = "yes" if "TANDEM" in arow[fidx["ALT"]] else ("no" if typ in ("DEL", "DUP") else "n/a")

        out_rows.append([
            gene, famA, famB, chrom, start, end, typ, gt, recip,
            arow[fidx["SV_length"]], arow[fidx["Location"]], tandem,
            arow[fidx["HI"]], arow[fidx["TS"]], arow[fidx["DDD_HI_percent"]],
            arow[fidx["GnomAD_pLI"]], arow[fidx["LOEUF_bin"]],
            arow[fidx["Overlapped_CDS_percent"]], arow[fidx["Frameshift"]],
            arow[fidx["Dist_nearest_SS"]], arow[fidx["Nearest_SS_type"]],
            srow[sidx["SpliceAI_max_score"]] if srow else "",
            srow[sidx["GREENDB_type"]] if srow else "",
            srow[sidx["GREENDB_constraint"]] if srow else "",
            arow[fidx["OMIM_phenotype"]],
            nfam, phA, phB, reason,
        ])

    wb = Workbook()
    ws = wb.active
    ws.title = "Genewide_overlap_candidates"
    ws.append(OUT_COLS)
    header_font = Font(bold=True)
    header_fill = PatternFill(start_color="D9E1F2", end_color="D9E1F2", fill_type="solid")
    for c in range(1, len(OUT_COLS) + 1):
        cell = ws.cell(row=1, column=c)
        cell.font = header_font
        cell.fill = header_fill
    ws.freeze_panes = "A2"
    for row in out_rows:
        ws.append(row)
    for c in range(1, len(OUT_COLS) + 1):
        col_letter = get_column_letter(c)
        width = max(len(str(OUT_COLS[c - 1])), *(len(str(row[c - 1])) for row in out_rows)) if out_rows else len(OUT_COLS[c - 1])
        ws.column_dimensions[col_letter].width = min(width + 2, 60)

    wb.save(args.out)
    print(f"Workbook written: {args.out} ({len(out_rows)} rows)", file=sys.stderr)


if __name__ == "__main__":
    main()
