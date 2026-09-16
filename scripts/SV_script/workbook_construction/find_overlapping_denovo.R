#!/usr/bin/env Rscript
#
# Finds de novo SVs recurring across probands: groups by gene, then checks
# which hits actually overlap in coordinates and share the same GT (same
# breakpoints but different zygosity isn't the same event). Writes a
# 2-sheet workbook: Gene_recurrence (genes hit in >=2 families) and
# Overlapping_SV_pairs (actual coordinate-overlapping pairs, split into
# same/different proband phenotype).
#
# Possible_artefact fires only when a pair is an EXACT coordinate match AND
# that same locus recurs in >=3 families genome-wide - independent de novo
# events landing on byte-identical breakpoints 3+ times is implausible, so
# this combination reads as a calling artefact/repeat hotspot rather than
# real recurrence. Both conditions are required (AND, not OR): 3+ family
# recurrence without exact agreement, or a 100% match in just 2 families,
# is left for manual review instead.
#
# Usage:
#   Rscript find_overlapping_denovo.R --out <report.xlsx> [--pheno <diagnosis.csv>] <family1_allgenes.tsv> <family2_allgenes.tsv> ...

library(optparse)
library(readr)
library(dplyr)
library(openxlsx)

get_script_dir = function() {
  file_arg = grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg) == 1) dirname(normalizePath(sub("^--file=", "", file_arg))) else getwd()
}
default_pheno = normalizePath(file.path(get_script_dir(), "..", "..", "..", "data",
                                         "DRAGEN_complete_trios_plus_siblings.csv"), mustWork = FALSE)

option_list = list(
  make_option("--out", type = "character"),
  make_option("--pheno", type = "character", default = default_pheno,
              help = "CSV with FamilyNo/ICONIK_proband_diagnosis/ICONIK_Family_member columns [default: %default]")
)
parsed = parse_args(OptionParser(option_list = option_list), positional_arguments = TRUE)
opt = parsed$options
tsv_paths = parsed$args

if (length(tsv_paths) == 0) {
  stop("No family TSV files given.")
}

# per-family proband phenotype (ICONIK clinical diagnosis: DCM/HCM/LVNC/...)
pheno_lookup = read_csv(opt$pheno, col_types = cols(.default = "c")) %>%
  filter(ICONIK_Family_member == "Proband") %>%
  mutate(Family = sub("^([0-9]+).*", "\\1", FamilyNo)) %>%
  group_by(Family) %>%
  summarise(Phenotype = dplyr::first(ICONIK_proband_diagnosis), .groups = "drop")

AF_COLS = c(DEL = "B_loss_AFmax", DUP = "B_gain_AFmax", INS = "B_ins_AFmax", INV = "B_inv_AFmax")

# AnnotSV's dosage-sensitivity/haploinsufficiency columns
HI_COLS = c("HI", "TS", "DDD_HI_percent", "GnomAD_pLI", "LOEUF_bin")

# min families sharing an exact locus to flag as a recurrent artefact (see header)
ARTEFACT_MIN_FAMILIES = 3

all_svs = bind_rows(lapply(tsv_paths, function(p) {
  tbl = read_tsv(p, col_types = cols(.default = "c"))
  # proband GT column name embeds the sample ID (Proband_<SM>_GT) - find it by pattern
  gt_col = grep("^Proband_.*_GT$", colnames(tbl), value = TRUE)
  if (length(gt_col) != 1) {
    stop(sprintf("Expected exactly one Proband_*_GT column in %s, found %d", p, length(gt_col)))
  }
  tbl$Proband_GT = tbl[[gt_col]]
  tbl[, c("Family", "Proband_Sex", "Proband_GT", "Gene_name", "SV_chrom", "SV_start", "SV_end", "SV_length", "SV_type",
          "B_loss_AFmax", "B_gain_AFmax", "B_ins_AFmax", "B_inv_AFmax", HI_COLS)]
}))
all_svs$SV_start = as.numeric(all_svs$SV_start)
all_svs$SV_end = as.numeric(all_svs$SV_end)
all_svs$SV_length = as.numeric(all_svs$SV_length)
all_svs = all_svs[!is.na(all_svs$SV_start) & !is.na(all_svs$SV_end), ]
# population AF for this SV's own type, whichever B_*_AFmax column applies
all_svs$AFmax = mapply(function(typ, row_idx) {
  col = AF_COLS[[typ]]
  if (is.null(col)) return(NA_character_)
  all_svs[[col]][row_idx]
}, all_svs$SV_type, seq_len(nrow(all_svs)))
# INS coordinates don't reflect true extent - use start+length as the effective
# end for overlap math instead of the near-meaningless raw SV_end
all_svs$SV_end_eff = ifelse(all_svs$SV_type == "INS" & !is.na(all_svs$SV_length),
                             all_svs$SV_start + abs(all_svs$SV_length) - 1,
                             all_svs$SV_end)

# any gene hit in >=2 different families
gene_recurrence = all_svs %>%
  group_by(Gene_name) %>%
  summarise(
    N_families = n_distinct(Family),
    N_SVs = n(),
    Families = paste(sort(unique(Family)), collapse = "; "),
    .groups = "drop"
  ) %>%
  filter(N_families >= 2) %>%
  arrange(desc(N_families))

message(sprintf("  %d genes hit in >=2 families (out of %d genes total)", nrow(gene_recurrence), n_distinct(all_svs$Gene_name)))

# distinct families sharing the exact same locus for a gene - feeds the artefact flag
locus_family_counts = all_svs %>%
  group_by(Gene_name, SV_chrom, SV_start, SV_end_eff) %>%
  summarise(N_families_exact_locus = n_distinct(Family), .groups = "drop")

# within each recurrent gene, compare every pair of SVs from different families
overlap_rows = list()

for (gene in gene_recurrence$Gene_name) {
  svs = all_svs[all_svs$Gene_name == gene, ]
  if (nrow(svs) < 2) next

  pairs = combn(nrow(svs), 2)
  for (k in seq_len(ncol(pairs))) {
    i = pairs[1, k]
    j = pairs[2, k]
    if (svs$Family[i] == svs$Family[j]) next  # cross-family pairs only
    if (svs$SV_chrom[i] != svs$SV_chrom[j]) next  # same gene symbol, different chrom (paralogue) - not a real overlap
    if (is.na(svs$Proband_GT[i]) || is.na(svs$Proband_GT[j]) || svs$Proband_GT[i] != svs$Proband_GT[j]) next  # different zygosity isn't the same event

    overlap_bp = min(svs$SV_end_eff[i], svs$SV_end_eff[j]) - max(svs$SV_start[i], svs$SV_start[j]) + 1
    if (overlap_bp <= 0) next

    len_i = abs(svs$SV_end_eff[i] - svs$SV_start[i]) + 1
    len_j = abs(svs$SV_end_eff[j] - svs$SV_start[j]) + 1
    pct_i = round(100 * overlap_bp / len_i, 1)
    pct_j = round(100 * overlap_bp / len_j, 1)
    reciprocal_pct = min(pct_i, pct_j)

    # dosage-sensitivity annotations should agree (same gene); fall back to whichever is non-NA
    hi_vals = setNames(lapply(HI_COLS, function(col) {
      vi = svs[[col]][i]; vj = svs[[col]][j]
      val = if (!is.na(vi)) vi else vj
      if (col %in% c("DDD_HI_percent", "GnomAD_pLI", "LOEUF_bin")) as.numeric(val) else val
    }), HI_COLS)

    afmax_i = as.numeric(svs$AFmax[i])
    afmax_j = as.numeric(svs$AFmax[j])

    # artefact flag - see header comment (exact match AND >=3-family recurrence)
    is_exact_match = (svs$SV_start[i] == svs$SV_start[j]) && (svs$SV_end_eff[i] == svs$SV_end_eff[j])
    n_families_exact_locus = if (is_exact_match) {
      locus_family_counts$N_families_exact_locus[
        locus_family_counts$Gene_name == gene &
        locus_family_counts$SV_chrom == svs$SV_chrom[i] &
        locus_family_counts$SV_start == svs$SV_start[i] &
        locus_family_counts$SV_end_eff == svs$SV_end_eff[i]]
    } else 0L

    possible_artefact = is_exact_match && n_families_exact_locus >= ARTEFACT_MIN_FAMILIES
    artefact_reason = if (possible_artefact) sprintf(
      "exact identical breakpoints (chr%s:%d-%d) shared by %d unrelated families for %s - probability of independent de novo recurrence at byte-identical coordinates is negligible; consistent with a confirmed recurrent locus-specific calling artefact or repeat-mediated hotspot",
      svs$SV_chrom[i], svs$SV_start[i], svs$SV_end_eff[i], n_families_exact_locus, gene) else NA_character_

    overlap_rows[[length(overlap_rows) + 1]] = data.frame(
      Gene_name = gene,
      Family_A = svs$Family[i], Sex_A = svs$Proband_Sex[i], Chrom_A = svs$SV_chrom[i], Start_A = svs$SV_start[i], End_A = svs$SV_end[i], Type_A = svs$SV_type[i], GT_A = svs$Proband_GT[i],
      Family_B = svs$Family[j], Sex_B = svs$Proband_Sex[j], Chrom_B = svs$SV_chrom[j], Start_B = svs$SV_start[j], End_B = svs$SV_end[j], Type_B = svs$SV_type[j], GT_B = svs$Proband_GT[j],
      Overlap_bp = overlap_bp,
      Pct_overlap_of_A = pct_i,
      Pct_overlap_of_B = pct_j,
      Reciprocal_pct_overlap = reciprocal_pct,
      AFmax_A = afmax_i, AFmax_B = afmax_j,
      HI = hi_vals[["HI"]], TS = hi_vals[["TS"]], DDD_HI_percent = hi_vals[["DDD_HI_percent"]],
      GnomAD_pLI = hi_vals[["GnomAD_pLI"]], LOEUF_bin = hi_vals[["LOEUF_bin"]],
      Possible_artefact = ifelse(possible_artefact, "Yes", "No"),
      Artefact_reason = artefact_reason,
      stringsAsFactors = FALSE
    )
  }
}

overlapping_pairs = if (length(overlap_rows) > 0) {
  bind_rows(overlap_rows) %>%
    left_join(pheno_lookup %>% rename(Family_A = Family, Phenotype_A = Phenotype), by = "Family_A") %>%
    left_join(pheno_lookup %>% rename(Family_B = Family, Phenotype_B = Phenotype), by = "Family_B") %>%
    mutate(Same_phenotype = !is.na(Phenotype_A) & !is.na(Phenotype_B) & Phenotype_A == Phenotype_B) %>%
    arrange(desc(Reciprocal_pct_overlap))
} else {
  data.frame(Gene_name = character(),
             Family_A = character(), Sex_A = character(), Chrom_A = character(), Start_A = numeric(), End_A = numeric(), Type_A = character(), GT_A = character(),
             Family_B = character(), Sex_B = character(), Chrom_B = character(), Start_B = numeric(), End_B = numeric(), Type_B = character(), GT_B = character(),
             Overlap_bp = numeric(), Pct_overlap_of_A = numeric(), Pct_overlap_of_B = numeric(),
             Reciprocal_pct_overlap = numeric(),
             AFmax_A = numeric(), AFmax_B = numeric(),
             HI = character(), TS = character(), DDD_HI_percent = numeric(),
             GnomAD_pLI = numeric(), LOEUF_bin = numeric(),
             Possible_artefact = character(), Artefact_reason = character(),
             Phenotype_A = character(), Phenotype_B = character(), Same_phenotype = logical())
}

message(sprintf("  %d overlapping SV pairs found across %d recurrent genes", nrow(overlapping_pairs), nrow(gene_recurrence)))

# AF for genes where every pair is a 100% reciprocal overlap (same SV shared
# identically) - partial-overlap genes left NA for manual inspection
full100_genes = overlapping_pairs %>%
  group_by(Gene_name) %>%
  summarise(all_100 = all(Reciprocal_pct_overlap == 100), .groups = "drop") %>%
  filter(all_100) %>%
  pull(Gene_name)

af_lookup = all_svs %>%
  filter(Gene_name %in% full100_genes) %>%
  group_by(Gene_name) %>%
  summarise(Allele_frequency = as.numeric(dplyr::first(AFmax)), .groups = "drop")

gene_recurrence = gene_recurrence %>%
  left_join(af_lookup, by = "Gene_name") %>%
  relocate(Allele_frequency, .after = Gene_name)

# Overlapping_SV_pairs keeps matching-diagnosis pairs; the rest move to Different_phenotype
same_phenotype_pairs = overlapping_pairs %>% filter(Same_phenotype)
diff_phenotype_pairs = overlapping_pairs %>% filter(!Same_phenotype)

message(sprintf("  %d/%d overlapping pairs have matching proband phenotype (rest -> Different_phenotype sheet)",
                 nrow(same_phenotype_pairs), nrow(overlapping_pairs)))

wb = createWorkbook()
header_style = createStyle(textDecoration = "bold", fgFill = "#D9E1F2", border = "Bottom")

write_sheet = function(tbl, sheet_name) {
  addWorksheet(wb, sheet_name)
  writeData(wb, sheet_name, tbl, headerStyle = header_style)
  freezePane(wb, sheet_name, firstRow = TRUE)
  if (nrow(tbl) > 0) setColWidths(wb, sheet_name, cols = seq_len(ncol(tbl)), widths = "auto")
}

write_sheet(gene_recurrence, "Gene_recurrence")
write_sheet(same_phenotype_pairs, "Overlapping_SV_pairs")
write_sheet(diff_phenotype_pairs, "Different_phenotype")

saveWorkbook(wb, opt$out, overwrite = TRUE)
message(sprintf("Workbook written: %s", opt$out))
