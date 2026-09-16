#!/usr/bin/env Rscript
#
# Combines proband + sibling per-family SV TSVs into one workbook: a
# "Proband" sheet, a "Sibling" sheet, and optionally an "Overlap" sheet of
# same-family proband-vs-sibling SVs sharing a gene with overlapping
# coordinates.
#
# Usage:
#   Rscript build_proband_sibling_workbook.R \
#       --proband <357_proband.tsv> <45_proband.tsv> ... \
#       --sibling <357_sibling.tsv> <45_sibling.tsv> ... \
#       --out results/workbooks/inherited_R135_proband_sibling.xlsx \
#       [--overlap]

library(optparse)
library(readr)
library(dplyr)
library(openxlsx)

option_list = list(
  make_option("--proband", type = "character", action = "store", dest = "proband_raw"),
  make_option("--sibling", type = "character", action = "store", dest = "sibling_raw"),
  make_option("--out", type = "character"),
  make_option("--overlap", action = "store_true", default = FALSE)
)

# optparse doesn't support multi-value options - parse --proband/--sibling/--out/--overlap manually
raw_args = commandArgs(trailingOnly = TRUE)

split_args = function(args) {
  groups = list(proband = character(0), sibling = character(0), out = character(0), overlap = FALSE)
  current = NULL
  for (a in args) {
    if (a == "--proband") { current = "proband"; next }
    if (a == "--sibling") { current = "sibling"; next }
    if (a == "--out")     { current = "out"; next }
    if (a == "--overlap") { groups$overlap = TRUE; current = NULL; next }
    if (is.null(current)) stop(sprintf("Unexpected argument '%s' before any --proband/--sibling/--out flag", a))
    groups[[current]] = c(groups[[current]], a)
  }
  groups
}

opt = split_args(raw_args)
if (length(opt$proband) == 0) stop("No --proband TSV files given.")
if (length(opt$sibling) == 0) stop("No --sibling TSV files given.")
if (length(opt$out) != 1) stop("Exactly one --out path required.")

AF_COLS = c(DEL = "B_loss_AFmax", DUP = "B_gain_AFmax", INS = "B_ins_AFmax", INV = "B_inv_AFmax")

read_role_tsvs = function(paths) {
  bind_rows(lapply(paths, function(p) {
    if (!file.exists(p)) {
      message(sprintf("  WARNING: %s not found, skipping", p))
      return(NULL)
    }
    read_tsv(p, col_types = cols(.default = "c"))
  }))
}

proband_tbl = read_role_tsvs(opt$proband)
sibling_tbl = read_role_tsvs(opt$sibling)

message(sprintf("Proband: %d variants across %d families", nrow(proband_tbl), n_distinct(proband_tbl$Family)))
message(sprintf("Sibling: %d variants across %d families", nrow(sibling_tbl), n_distinct(sibling_tbl$Family)))

# both TSVs write the child slot as Proband_* - rename on the sibling copy for clarity
sibling_tbl_display = sibling_tbl
colnames(sibling_tbl_display) = sub("^Proband_", "Sibling_", colnames(sibling_tbl_display))

wb = createWorkbook()
header_style = createStyle(textDecoration = "bold", fgFill = "#D9E1F2", border = "Bottom")

write_sheet = function(tbl, sheet_name) {
  addWorksheet(wb, sheet_name)
  writeData(wb, sheet_name, tbl, headerStyle = header_style)
  freezePane(wb, sheet_name, firstRow = TRUE)
  if (nrow(tbl) > 0) setColWidths(wb, sheet_name, cols = seq_len(ncol(tbl)), widths = "auto")
}

write_sheet(proband_tbl, "Proband")
write_sheet(sibling_tbl_display, "Sibling")

if (opt$overlap) {

  prep = function(tbl) {
    tbl = tbl[, intersect(colnames(tbl), c("Family", "Gene_name", "SV_chrom", "SV_start", "SV_end",
                                            "SV_type", AF_COLS, "HI", "TS", "DDD_HI_percent",
                                            "GnomAD_pLI", "LOEUF_bin"))]
    tbl$SV_start = as.numeric(tbl$SV_start)
    tbl$SV_end = as.numeric(tbl$SV_end)
    tbl[!is.na(tbl$SV_start) & !is.na(tbl$SV_end), ]
  }

  p = prep(proband_tbl)
  s = prep(sibling_tbl)

  overlap_rows = list()
  common_families = intersect(unique(p$Family), unique(s$Family))

  for (fam in common_families) {
    p_fam = p[p$Family == fam, ]
    s_fam = s[s$Family == fam, ]
    common_genes = intersect(unique(p_fam$Gene_name), unique(s_fam$Gene_name))

    for (gene in common_genes) {
      p_svs = p_fam[p_fam$Gene_name == gene, ]
      s_svs = s_fam[s_fam$Gene_name == gene, ]

      for (i in seq_len(nrow(p_svs))) {
        for (j in seq_len(nrow(s_svs))) {
          if (p_svs$SV_chrom[i] != s_svs$SV_chrom[j]) next

          overlap_bp = min(p_svs$SV_end[i], s_svs$SV_end[j]) - max(p_svs$SV_start[i], s_svs$SV_start[j]) + 1
          if (overlap_bp <= 0) next

          len_p = abs(p_svs$SV_end[i] - p_svs$SV_start[i]) + 1
          len_s = abs(s_svs$SV_end[j] - s_svs$SV_start[j]) + 1
          pct_p = round(100 * overlap_bp / len_p, 1)
          pct_s = round(100 * overlap_bp / len_s, 1)

          overlap_rows[[length(overlap_rows) + 1]] = data.frame(
            Family = fam, Gene_name = gene,
            Proband_Chrom = p_svs$SV_chrom[i], Proband_Start = p_svs$SV_start[i],
            Proband_End = p_svs$SV_end[i], Proband_Type = p_svs$SV_type[i],
            Sibling_Chrom = s_svs$SV_chrom[j], Sibling_Start = s_svs$SV_start[j],
            Sibling_End = s_svs$SV_end[j], Sibling_Type = s_svs$SV_type[j],
            Overlap_bp = overlap_bp,
            Pct_overlap_of_proband = pct_p,
            Pct_overlap_of_sibling = pct_s,
            Reciprocal_pct_overlap = min(pct_p, pct_s),
            stringsAsFactors = FALSE
          )
        }
      }
    }
  }

  overlap_tbl = if (length(overlap_rows) > 0) {
    bind_rows(overlap_rows) %>% arrange(desc(Reciprocal_pct_overlap))
  } else {
    data.frame(Family = character(), Gene_name = character(),
               Proband_Chrom = character(), Proband_Start = numeric(), Proband_End = numeric(), Proband_Type = character(),
               Sibling_Chrom = character(), Sibling_Start = numeric(), Sibling_End = numeric(), Sibling_Type = character(),
               Overlap_bp = numeric(), Pct_overlap_of_proband = numeric(), Pct_overlap_of_sibling = numeric(),
               Reciprocal_pct_overlap = numeric())
  }

  message(sprintf("Overlap: %d proband/sibling SV pairs sharing a gene with overlapping coordinates", nrow(overlap_tbl)))
  write_sheet(overlap_tbl, "Overlap")
}

saveWorkbook(wb, opt$out, overwrite = TRUE)
message(sprintf("Workbook written: %s", opt$out))

# openxlsx 4.2.4 can leave dangling drawing/vmlDrawing rels that make Excel
# refuse to open the file. If that happens again, strip them with:
#   xml <- readLines(rels_file); gsub('<Relationship [^>]*Type="[^"]*/(drawing|vmlDrawing)"[^>]*/>', "", xml)
# (and the matching <Override PartName="/xl/drawings/..."> in [Content_Types].xml) then re-zip.
