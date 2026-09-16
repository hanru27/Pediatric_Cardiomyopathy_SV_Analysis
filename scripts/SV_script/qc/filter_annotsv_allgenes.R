#!/usr/bin/env Rscript
#
# Same as filter_annotsv_panel.R but keeps every gene (genome-wide candidate
# list instead of the R135 panel). Attaches per-sample genotypes and drops
# proband calls where FT != PASS.
#
# Usage:
#   Rscript filter_annotsv_allgenes.R \
#       --annotsv    <denovo_annotated.tsv> \
#       --gt         <denovo_genotypes.tsv> \
#       --family     <FAMILY_ID> \
#       --proband    <PROBAND_SM> \
#       --pedigree   <iconik_trioPedigree.csv> \
#       --proband_id <PROBAND_ID (trios.txt)> \
#       --mother_id  <MOTHER_ID (trios.txt)> \
#       --father_id  <FATHER_ID (trios.txt)> \
#       --out        <output.tsv>

library(optparse)
library(readr)
library(dplyr)

option_list = list(
  make_option("--annotsv",    type = "character"),
  make_option("--gt",         type = "character"),
  make_option("--family",     type = "character"),
  make_option("--proband",    type = "character"),
  make_option("--pedigree",   type = "character"),
  make_option("--proband_id", type = "character"),
  make_option("--mother_id",  type = "character"),
  make_option("--father_id",  type = "character"),
  make_option("--out",        type = "character")
)
opt = parse_args(OptionParser(option_list = option_list))

strip_chr = function(chrom) {
  sub("^chr", "", chrom, ignore.case = TRUE)
}

# trios.txt IDs sometimes have a trailing _D<n> the pedigree file doesn't so strip before matching
strip_id_suffix = function(id) {
  sub("_D[0-9]+$", "", id)
}

# fallback key for pedigree rows stored as a bare sample ID instead of "Ware_<id>_01[_D<n>]"
core_id = function(id) {
  id = strip_id_suffix(id)
  id = sub("^Ware_", "", id)
  sub("_[0-9]+$", "", id)
}

# pedigree lookup: proband sex, and affected status for all three
pedigree_tbl = read_csv(opt$pedigree, show_col_types = FALSE)

lookup_pedigree = function(individual_id, field) {
  row = pedigree_tbl[pedigree_tbl$`Individual ID` == strip_id_suffix(individual_id), ]
  if (nrow(row) == 0) {
    row = pedigree_tbl[pedigree_tbl$`Individual ID` == core_id(individual_id), ]
  }
  if (nrow(row) == 0) return("Unknown")
  row[[field]][1]
}

proband_sex      = lookup_pedigree(opt$proband_id, "Sex")
proband_affected = lookup_pedigree(opt$proband_id, "Affected Status")
mother_affected   = lookup_pedigree(opt$mother_id,  "Affected Status")
father_affected   = lookup_pedigree(opt$father_id,  "Affected Status")

# genotype table (proband, mother, father)
gt_tbl = read_tsv(opt$gt, col_types = cols(.default = "c"))
gt_tbl$join_chrom = strip_chr(gt_tbl$CHROM)
gt_tbl$join_key = paste(gt_tbl$join_chrom, gt_tbl$POS, gt_tbl$SVTYPE, sep = "|")

gt_extra_cols = setdiff(colnames(gt_tbl), c("CHROM", "POS", "SVTYPE", "join_chrom", "join_key"))
proband_ft_col = paste0("Proband_", opt$proband, "_FT")
if (!proband_ft_col %in% gt_extra_cols) {
  stop(sprintf("Column '%s' not found in %s — check --proband matches the sample name used in the genotype TSV header", proband_ft_col, opt$gt))
}

annotsv_tbl = read_tsv(opt$annotsv, col_types = cols(.default = "c"))

# keep only "split" rows (one per SV per gene) - "full" rows (all overlapped genes combined
# into one Gene_name) would otherwise duplicate an SV that overlaps more than one gene
n_before_split_only = nrow(annotsv_tbl)
annotsv_tbl = annotsv_tbl[annotsv_tbl$Annotation_mode == "split", ]
message(sprintf("  split-mode only (dropped full-mode summary rows): %d -> %d rows", n_before_split_only, nrow(annotsv_tbl)))

annotsv_tbl$join_chrom = strip_chr(annotsv_tbl$SV_chrom)
annotsv_tbl$join_key = paste(annotsv_tbl$join_chrom, annotsv_tbl$SV_start, annotsv_tbl$SV_type, sep = "|")

gt_by_key = gt_tbl[match(annotsv_tbl$join_key, gt_tbl$join_key), gt_extra_cols]
colnames(gt_by_key) = gt_extra_cols
gt_by_key[is.na(gt_by_key)] = "."
annotsv_tbl = bind_cols(annotsv_tbl, gt_by_key)

# proband de novo confidence: sample-level filter must be PASS
n_before_ft = nrow(annotsv_tbl)
annotsv_tbl = annotsv_tbl[annotsv_tbl[[proband_ft_col]] == "PASS", ]
message(sprintf("  proband FT=PASS filter: %d -> %d rows", n_before_ft, nrow(annotsv_tbl)))

# no gene panel filter here - every gene is kept
annotsv_tbl$Family = opt$family
annotsv_tbl$Proband_Sex = proband_sex
annotsv_tbl$Proband_Affected = proband_affected
annotsv_tbl$Mother_Affected = mother_affected
annotsv_tbl$Father_Affected = father_affected

# dedup: AnnotSV writes one row per overlapped transcript, so a multi-transcript
# gene can repeat the same SV - collapse rows with identical SV + gene
n_before_dedup = nrow(annotsv_tbl)
dedup_key = paste(annotsv_tbl$Family, annotsv_tbl$SV_start, annotsv_tbl$SV_end,
                   annotsv_tbl$SV_length, annotsv_tbl$SV_type, annotsv_tbl$Gene_name, sep = "|")
annotsv_tbl = annotsv_tbl[!duplicated(dedup_key), ]
message(sprintf("  duplicate calls removed: %d -> %d rows", n_before_dedup, nrow(annotsv_tbl)))

# rank by AnnotSV score, descending
rank_score = suppressWarnings(as.numeric(annotsv_tbl$AnnotSV_ranking_score))
rank_score[is.na(rank_score)] = 0
annotsv_tbl = annotsv_tbl[order(-rank_score), ]

# column order: Family/sex/affected, genotypes, then AnnotSV columns (incl. Gene_name)
leading = c("Family", "Proband_Sex", "Proband_Affected", "Mother_Affected", "Father_Affected", gt_extra_cols)
rest = setdiff(colnames(annotsv_tbl), c(leading, "join_chrom", "join_key"))
annotsv_tbl = annotsv_tbl[, c(leading, rest)]

write_tsv(annotsv_tbl, opt$out)
message(sprintf("  %d SVs (all genes) -> %s", nrow(annotsv_tbl), opt$out))
