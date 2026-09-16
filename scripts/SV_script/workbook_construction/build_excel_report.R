#!/usr/bin/env Rscript
#
# Combine per-family R135 green-gene de novo SV tables into a single Excel
# workbook, one sheet per family/proband (rather than one flat combined
# sheet), so each trio can be reviewed on its own tab.
#
# Usage:
#   Rscript build_excel_report.R --out <report.xlsx> <family1.tsv> <family2.tsv> ...

library(optparse)
library(readr)
library(openxlsx)

option_list = list(
  make_option("--out", type = "character")
)
parsed = parse_args(OptionParser(option_list = option_list), positional_arguments = TRUE)
opt = parsed$options
tsv_paths = parsed$args

if (length(tsv_paths) == 0) {
  stop("No family TSV files given.")
}

wb = createWorkbook()
header_style = createStyle(textDecoration = "bold", fgFill = "#D9E1F2", border = "Bottom")
used_sheet_names = character(0)

for (tsv_path in tsv_paths) {
  if (!file.exists(tsv_path)) {
    message(sprintf("  WARNING: %s not found, skipping", tsv_path))
    next
  }

  tbl = read_tsv(tsv_path, col_types = cols(.default = "c"))

  family_id = if (nrow(tbl) > 0) tbl$Family[1] else basename(dirname(dirname(dirname(tsv_path))))
  sheet_name = paste0("Family_", family_id)
  sheet_name = substr(sheet_name, 1, 31)
  if (sheet_name %in% used_sheet_names) {
    sheet_name = substr(paste0(sheet_name, "_", length(used_sheet_names)), 1, 31)
  }
  used_sheet_names = c(used_sheet_names, sheet_name)

  addWorksheet(wb, sheet_name)
  writeData(wb, sheet_name, tbl, headerStyle = header_style)
  freezePane(wb, sheet_name, firstRow = TRUE)
  if (nrow(tbl) > 0) {
    setColWidths(wb, sheet_name, cols = seq_len(ncol(tbl)), widths = "auto")
  }

  message(sprintf("  %s: %d variants -> sheet '%s'", basename(tsv_path), nrow(tbl), sheet_name))
}

saveWorkbook(wb, opt$out, overwrite = TRUE)
message(sprintf("Workbook written: %s (%d sheets)", opt$out, length(used_sheet_names)))
