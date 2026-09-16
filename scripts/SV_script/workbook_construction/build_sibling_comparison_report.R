#!/usr/bin/env Rscript
#
# For each family, combines the proband's and sibling's AnnotSV-filtered SV
# table into one workbook: a proband sheet, a sibling sheet, and (with
# --overlap) a sheet of proband/sibling SV pairs overlapping in coordinates
# (same chrom + SV_type), with reciprocal % overlap.
#
# INS breakpoints are a single reference point, not a span, so coordinate
# overlap is meaningless there - INS pairs get pairwise ALT sequence
# identity (edit-distance based) instead.
#
# Usage:
#   Rscript build_sibling_comparison_report.R --out results/workbooks/sibling_comparison.xlsx [--overlap] \
#       <family1> <proband1.tsv> <sibling1.tsv> [<family2> <proband2.tsv> <sibling2.tsv> ...]

library(optparse)
library(readr)
library(openxlsx)

option_list = list(
  make_option("--out", type = "character"),
  make_option("--overlap", action = "store_true", default = FALSE,
              help = "Add a per-family Overlap sheet comparing proband vs sibling SVs")
)
parsed = parse_args(OptionParser(option_list = option_list), positional_arguments = TRUE)
opt = parsed$options
args = parsed$args

if (length(args) == 0 || length(args) %% 3 != 0) {
  stop("Expected triples of: <family> <proband_tsv> <sibling_tsv>")
}

read_sv_tbl = function(path) {
  if (!file.exists(path)) {
    message(sprintf("  WARNING: %s not found, using empty table", path))
    return(NULL)
  }
  read_tsv(path, col_types = cols(.default = "c"))
}

empty_overlap = data.frame(
  Gene_name = character(), Chrom = character(),
  Proband_Start = numeric(), Proband_End = numeric(), Proband_Type = character(),
  Sibling_Start = numeric(), Sibling_End = numeric(), Sibling_Type = character(),
  Overlap_bp = numeric(), Pct_overlap_of_proband = numeric(),
  Pct_overlap_of_sibling = numeric(), Reciprocal_pct_overlap = numeric(),
  stringsAsFactors = FALSE
)

# edit-distance-based identity between two ALT sequences, for INS pairs
seq_identity_pct = function(seq1, seq2) {
  if (is.na(seq1) || is.na(seq2) || nchar(seq1) == 0 || nchar(seq2) == 0) return(NA_real_)
  edit_dist = adist(seq1, seq2)[1, 1]
  round(100 * (1 - edit_dist / max(nchar(seq1), nchar(seq2))), 1)
}

# every proband SV vs every sibling SV in the family, same chrom + SV_type required
compute_overlap = function(proband_tbl, sibling_tbl) {
  if (is.null(proband_tbl) || is.null(sibling_tbl) || nrow(proband_tbl) == 0 || nrow(sibling_tbl) == 0) {
    return(empty_overlap)
  }
  p = proband_tbl
  s = sibling_tbl
  p$SV_start = as.numeric(p$SV_start); p$SV_end = as.numeric(p$SV_end)
  s$SV_start = as.numeric(s$SV_start); s$SV_end = as.numeric(s$SV_end)

  rows = list()
  for (i in seq_len(nrow(p))) {
    cand = which(s$SV_chrom == p$SV_chrom[i] & s$SV_type == p$SV_type[i])
    for (j in cand) {
      overlap_bp = min(p$SV_end[i], s$SV_end[j]) - max(p$SV_start[i], s$SV_start[j]) + 1
      if (is.na(overlap_bp) || overlap_bp <= 0) next

      if (p$SV_type[i] == "INS") {
        pct_p = seq_identity_pct(p$ALT[i], s$ALT[j])
        pct_s = pct_p
      } else {
        len_p = abs(p$SV_end[i] - p$SV_start[i]) + 1
        len_s = abs(s$SV_end[j] - s$SV_start[j]) + 1
        pct_p = round(100 * overlap_bp / len_p, 1)
        pct_s = round(100 * overlap_bp / len_s, 1)
      }

      rows[[length(rows) + 1]] = data.frame(
        Gene_name = p$Gene_name[i],
        Chrom = p$SV_chrom[i],
        Proband_Start = p$SV_start[i], Proband_End = p$SV_end[i], Proband_Type = p$SV_type[i],
        Sibling_Start = s$SV_start[j], Sibling_End = s$SV_end[j], Sibling_Type = s$SV_type[j],
        Overlap_bp = overlap_bp,
        Pct_overlap_of_proband = pct_p,
        Pct_overlap_of_sibling = pct_s,
        Reciprocal_pct_overlap = min(pct_p, pct_s),
        stringsAsFactors = FALSE
      )
    }
  }
  if (length(rows) == 0) return(empty_overlap)
  out = do.call(rbind, rows)
  out[order(-out$Reciprocal_pct_overlap), ]
}

wb = createWorkbook()
header_style = createStyle(textDecoration = "bold", fgFill = "#D9E1F2", border = "Bottom")

write_sheet = function(tbl, sheet_name) {
  sheet_name = substr(sheet_name, 1, 31)
  addWorksheet(wb, sheet_name)
  if (!is.null(tbl)) {
    writeData(wb, sheet_name, tbl, headerStyle = header_style)
    freezePane(wb, sheet_name, firstRow = TRUE)
    if (nrow(tbl) > 0) setColWidths(wb, sheet_name, cols = seq_len(ncol(tbl)), widths = "auto")
  }
  sheet_name
}

n_families = length(args) / 3
for (k in seq_len(n_families)) {
  family = args[(k - 1) * 3 + 1]
  proband_path = args[(k - 1) * 3 + 2]
  sibling_path = args[(k - 1) * 3 + 3]

  proband_tbl = read_sv_tbl(proband_path)
  sibling_tbl = read_sv_tbl(sibling_path)

  write_sheet(proband_tbl, sprintf("%s_Proband", family))
  write_sheet(sibling_tbl, sprintf("%s_Sibling", family))
  message(sprintf("  Family %s: proband %d variants, sibling %d variants",
                   family, if (is.null(proband_tbl)) 0 else nrow(proband_tbl),
                   if (is.null(sibling_tbl)) 0 else nrow(sibling_tbl)))

  if (opt$overlap) {
    overlap_tbl = compute_overlap(proband_tbl, sibling_tbl)
    write_sheet(overlap_tbl, sprintf("%s_Overlap", family))
    message(sprintf("  Family %s: %d overlapping proband/sibling SV pairs", family, nrow(overlap_tbl)))
  }
}

saveWorkbook(wb, opt$out, overwrite = TRUE)
message(sprintf("Workbook written: %s (%d sheets)", opt$out, length(wb$sheet_names)))
