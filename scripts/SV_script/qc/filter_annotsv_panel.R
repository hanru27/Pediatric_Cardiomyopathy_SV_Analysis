#!/usr/bin/env Rscript
#
# Filters AnnotSV output to R135 green-gene panel SVs, attaches per-sample
# genotypes (GT/GQ/FT/PR/SR), and joins in SpliceAI/GREEN-DB scores if given.
# Also drops proband calls where FT != PASS - tier-1's DQ=60 only catches
# Mendelian-inconsistent genotypes, not a bad proband-level filter call.
# GenCC columns are dropped (found unreliable for SV interpretation).
#
# Usage:
#   Rscript filter_annotsv_panel.R \
#       --annotsv    <denovo_annotated.tsv> \
#       --genes      <R135_green_genes.tsv> \
#       --gt         <denovo_genotypes.tsv> \
#       --family     <FAMILY_ID> \
#       --proband    <PROBAND_SM> \
#       --pedigree   <iconik_trioPedigree.csv> \
#       --proband_id <PROBAND_ID (trios.txt)> \
#       --mother_id  <MOTHER_ID (trios.txt)> \
#       --father_id  <FATHER_ID (trios.txt)> \
#       --splice_greendb <family_splice_greendb.tsv, optional> \
#       --out        <output.tsv>

library(optparse)
library(readr)
library(dplyr)
library(stringr)

option_list = list(
  make_option("--annotsv",    type = "character"),
  make_option("--genes",      type = "character"),
  make_option("--gt",         type = "character"),
  make_option("--family",     type = "character"),
  make_option("--proband",    type = "character"),
  make_option("--pedigree",   type = "character"),
  make_option("--proband_id", type = "character"),
  make_option("--mother_id",  type = "character"),
  make_option("--father_id",  type = "character"),
  make_option("--splice_greendb", type = "character", default = NA_character_),
  make_option("--out",        type = "character")
)
opt = parse_args(OptionParser(option_list = option_list))

GENCC_COLS = c("GenCC_disease", "GenCC_moi", "GenCC_classification", "GenCC_pmid")
SPLICE_GREENDB_COLS = c("SpliceAI_max_score", "SpliceAI_gene", "SpliceAI_effect", "SpliceAI_breakpoint",
                        "GREENDB_id", "GREENDB_type", "GREENDB_source", "GREENDB_genes", "GREENDB_constraint")

strip_chr = function(chrom) {
  sub("^chr", "", chrom, ignore.case = TRUE)
}

# trios.txt IDs sometimes have a trailing _D<n> the pedigree file doesn't - strip before matching
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

# green gene panel
green_genes = read_tsv(opt$genes, show_col_types = FALSE) %>%
  pull(Gene) %>%
  str_trim() %>%
  unique()

split_genes = function(gene_name) {
  str_trim(str_split(gene_name, "[/;,]")[[1]])
}

`%||%` = function(a, b) if (is.null(a) || is.na(a)) b else a

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
annotsv_tbl = annotsv_tbl[, !colnames(annotsv_tbl) %in% GENCC_COLS]

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

# SpliceAI + GREEN-DB, matched on the raw Manta SV ID
if (!is.na(opt$splice_greendb)) {
  splice_greendb_tbl = read_tsv(opt$splice_greendb, col_types = cols(.default = "c"))
  sg_by_id = splice_greendb_tbl[match(annotsv_tbl$ID, splice_greendb_tbl$SV_ID), SPLICE_GREENDB_COLS]
  colnames(sg_by_id) = SPLICE_GREENDB_COLS
  annotsv_tbl = bind_cols(annotsv_tbl, sg_by_id)
} else {
  annotsv_tbl[, SPLICE_GREENDB_COLS] = NA_character_
}

# proband de novo confidence: sample-level filter must be PASS
n_before_ft = nrow(annotsv_tbl)
annotsv_tbl = annotsv_tbl[annotsv_tbl[[proband_ft_col]] == "PASS", ]
message(sprintf("  proband FT=PASS filter: %d -> %d rows", n_before_ft, nrow(annotsv_tbl)))

# green gene panel filter
gene_hits = lapply(annotsv_tbl$Gene_name, function(gene_name) {
  parts = split_genes(gene_name %||% "")
  parts[parts %in% green_genes]
})

keep = lengths(gene_hits) > 0
annotsv_tbl = annotsv_tbl[keep, ]
annotsv_tbl$R135_green_genes = vapply(gene_hits[keep], paste, character(1), collapse = "; ")
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

# column order: Family/sex/affected, R135_green_genes, genotypes, SpliceAI/GREEN-DB, then AnnotSV columns
leading = c("Family", "Proband_Sex", "Proband_Affected", "Mother_Affected", "Father_Affected",
            "R135_green_genes", gt_extra_cols, SPLICE_GREENDB_COLS)
rest = setdiff(colnames(annotsv_tbl), c(leading, "join_chrom", "join_key"))
annotsv_tbl = annotsv_tbl[, c(leading, rest)]

write_tsv(annotsv_tbl, opt$out)
message(sprintf("  %d R135 green-gene SVs -> %s", nrow(annotsv_tbl), opt$out))
