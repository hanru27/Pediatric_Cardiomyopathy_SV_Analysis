#!/usr/bin/env Rscript
#
# Drops SVs above a benign population-AF threshold, checked only against the
# AF column matching the row's own SV_type (DEL -> B_loss_AFmax, DUP ->
# B_gain_AFmax, INS -> B_ins_AFmax) - a benign deletion's frequency says
# nothing about a duplication at the same locus. Inversions are out of scope
# (no B_inv_AFmax mapping); missing AF values pass, a value at the threshold
# exactly also passes.
#
# Usage:
#   Rscript filter_by_benign_af.R --threshold 0.01 --in <in.tsv> --out <out.tsv>

library(optparse)
library(readr)

option_list = list(
  make_option("--threshold", type = "double", default = 0.01),
  make_option("--in",  type = "character", dest = "in_path"),
  make_option("--out", type = "character")
)
opt = parse_args(OptionParser(option_list = option_list))

# B_inv_AFmax deliberately absent - inversions aren't analysed here
TYPE_AF_COL = c(DEL = "B_loss_AFmax", DUP = "B_gain_AFmax", INS = "B_ins_AFmax")

tbl = read_tsv(opt$in_path, col_types = cols(.default = "c"))

if (nrow(tbl) == 0) {
  write_tsv(tbl, opt$out)
  message(sprintf("  benign AF filter (threshold %.2f): 0 -> 0 rows", opt$threshold))
  quit(save = "no", status = 0)
}

af_col_per_row = TYPE_AF_COL[tbl$SV_type]
af_value = mapply(function(col, i) {
  if (is.na(col)) return(NA_real_)
  suppressWarnings(as.numeric(tbl[[col]][i]))
}, af_col_per_row, seq_len(nrow(tbl)))
fails = !is.na(af_value) & af_value > opt$threshold

n_before = nrow(tbl)
tbl = tbl[!fails, ]
message(sprintf("  benign AF filter (threshold %.2f, type-matched AF column, inversions excluded): %d -> %d rows", opt$threshold, n_before, nrow(tbl)))

write_tsv(tbl, opt$out)
