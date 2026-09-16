#!/usr/bin/env Rscript
# Per-sample SV size violin+points plots: raw vs tier1-filtered, and de novo
# vs inherited. Recreates the 2 combined PNGs in figure_analysis_2/SV/frequency_charts/.
# Requires: bcftools on PATH, R pkgs data.table/ggplot2/scales/patchwork.

library(data.table)
library(ggplot2)
library(scales)
library(patchwork)

base_dir = "/rds/general/user/rg3825/home/imperial_projects/summer_project/iconik_sv_analysis"
family_list = file.path(base_dir, "data/trios_all37.txt")
variants_dir = file.path(base_dir, "results/De_novo_calls")
out_dir = file.path(base_dir, "figure_analysis_2/SV/frequency_charts")
cache_dir = file.path(out_dir, "cache")
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

POINTS_PER_GROUP = 3000
VIOLIN_ADJUST = 1.8
REF_LINES = c(25000, 50)
sv_levels = c("DEL", "DUP", "INS")
families = as.character(fread(family_list, header = FALSE)$V1)

# one row per (site, sample-carrier) - a variant carried by 2 of 3 trio members counts twice
extract_persample = function(vcf_path, family) {
  raw = system2("bcftools", c("query", "-f", shQuote("%SVTYPE\t%INFO/SVLEN[\t%SAMPLE=%GT]\n"), shQuote(vcf_path)),
                stdout = TRUE, stderr = FALSE)
  if (!length(raw)) return(NULL)
  rows = lapply(strsplit(raw, "\t", fixed = TRUE), function(x) {
    svtype = x[1]
    svlen = x[2]
    gts = x[-(1:2)]
    gt_only = sub("^[^=]+=", "", gts)
    carriers = gts[grepl("1", gt_only) & !grepl("\\.", gt_only)]
    if (!length(carriers)) return(NULL)
    data.table(sample = paste0(family, ":", sub("=.*$", "", carriers)), svtype = svtype, svlen = svlen)
  })
  rbindlist(rows)
}

# pools one VCF tier across all 37 families, cached since bcftools extraction is slow
build_dt = function(vcf_name, tag) {
  cache_file = file.path(cache_dir, paste0("sv_persample_", tag, ".rds"))
  if (file.exists(cache_file)) return(readRDS(cache_file))
  parts = list()
  for (fam in families) {
    candidates = c(file.path(variants_dir, fam, "results/variants", vcf_name),
                   file.path(variants_dir, fam, "proband/results/variants", vcf_name))
    vcf_path = candidates[file.exists(candidates)][1]
    if (is.na(vcf_path)) {
      message("WARNING: no ", vcf_name, " for family ", fam, ", skipping")
      next
    }
    parts[[fam]] = extract_persample(vcf_path, fam)
  }
  dt = rbindlist(parts)
  dt[, svlen := abs(suppressWarnings(as.numeric(svlen)))]
  dt[svtype == "", svtype := "NA"]
  dt[, svtype := factor(svtype)]
  saveRDS(dt, cache_file)
  dt
}

clean_size_dt = function(dt) {
  size_dt = dt[!svtype %in% c("BND", "TRA") & is.finite(svlen) & svlen > 0]
  size_dt[, svtype := factor(as.character(svtype), levels = sv_levels)]
  size_dt[!is.na(svtype)]
}

# true log-position minor ticks (2,3,4...9x10^n per decade), all drawn at the same length
log_minor_breaks = function(limits) {
  vals = as.numeric(1:9 %o% 10^(0:ceiling(log10(limits[2]))))
  vals[vals >= limits[1] & vals <= limits[2]]
}

# top-right "Median Size & Total Count" inset box
median_box_label = function(size_dt) {
  stats = size_dt[, .(median_bp = median(svlen), n = .N), by = svtype][match(sv_levels, svtype)]
  header = sprintf("%-4s%12s%14s", "", "Median Size", "Total Count")
  rows = sprintf("%-4s%12s%14s", stats$svtype, format(round(stats$median_bp), big.mark = ","), format(stats$n, big.mark = ","))
  paste0(header, "\n", paste(rows, collapse = "\n"))
}

# violin + jittered subsample of points (full data still backs the violin/boxplot shape)
make_violin_points_plot = function(size_dt, title, y_limits, y_breaks, y_minor_breaks) {
  set.seed(1)
  pts = size_dt[, .SD[sample(.N, min(.N, POINTS_PER_GROUP))], by = svtype]
  ggplot(size_dt, aes(svtype, svlen, fill = svtype)) +
    geom_violin(scale = "width", trim = TRUE, adjust = VIOLIN_ADJUST, colour = "#222222") +
    geom_jitter(data = pts, width = 0.12, height = 0, size = 0.6, alpha = 0.35, colour = "black", shape = 16, show.legend = FALSE) +
    geom_hline(yintercept = REF_LINES, colour = "red", linetype = "dashed", linewidth = 0.4) +
    annotate("label", x = Inf, y = Inf, label = median_box_label(size_dt), hjust = 1, vjust = 1.05,
             size = 3.6, family = "mono", fill = "white", colour = "#222222") +
    scale_y_log10(labels = comma, breaks = y_breaks, minor_breaks = y_minor_breaks, limits = y_limits,
                  guide = guide_axis(minor.ticks = TRUE)) +
    labs(x = "SV type", y = "SV length (bp, log10)", title = title) +
    theme_bw() +
    theme(legend.position = "none", panel.grid.minor = element_blank(),
          axis.minor.ticks.length = rel(0.5), plot.title = element_text(hjust = 0.5))
}

# shared y-axis (99.5th pctile of pooled raw+tier1, snapped to a clean log decade) for every plot
dt_raw = build_dt("diploidSV.vcf.gz", "raw")
dt_tier1 = build_dt("diploidSV.tier1.vcf.gz", "tier1")
size_raw = clean_size_dt(dt_raw)
size_tier1 = clean_size_dt(dt_tier1)

p995 = as.numeric(quantile(c(size_raw$svlen, size_tier1$svlen), 0.995))
candidate_tops = c(1e4, 3e4, 1e5, 3e5, 1e6, 3e6, 1e7)
y_top = min(candidate_tops[candidate_tops >= p995])
y_breaks = sort(unique(c(10^(1:ceiling(log10(y_top))), REF_LINES)))
y_limits = c(50, y_top)
y_minor_breaks = log_minor_breaks(y_limits)

title_raw = "Violin plot of distribution of SV's across ICONIK cohort,\npre QC and filtering"
title_tier1 = "Violin plot of distribution of SV's across ICONIK cohort,\npost QC and filtering (Size, PASS, BND)"
p_raw = make_violin_points_plot(size_raw, title_raw, y_limits, y_breaks, y_minor_breaks)
p_tier1 = make_violin_points_plot(size_tier1, title_tier1, y_limits, y_breaks, y_minor_breaks)
ggsave(file.path(out_dir, "sv_size_violin_points_defaultcolors_combined.png"), p_raw | p_tier1, width = 15, height = 7, dpi = 200)

# de novo vs inherited, same fixed y-axis so it's comparable to the plot above
dt_denovo = build_dt("diploidSV.tier1_denovo_candidates.vcf.gz", "denovo")
dt_inherited = build_dt("diploidSV.tier1_inherited_candidates.vcf.gz", "inherited")
size_denovo = clean_size_dt(dt_denovo)
size_inherited = clean_size_dt(dt_inherited)

title_denovo = "Violin plot of distribution of SV's across ICONIK cohort\npost QC and filtering, de novo candidates only"
title_inherited = "Violin plot of distribution of SV's across ICONIK cohort\npost QC and filtering, inherited candidates only"
p_denovo = make_violin_points_plot(size_denovo, title_denovo, y_limits, y_breaks, y_minor_breaks)
p_inherited = make_violin_points_plot(size_inherited, title_inherited, y_limits, y_breaks, y_minor_breaks)
ggsave(file.path(out_dir, "sv_size_violin_points_defaultcolors_denovo_vs_inherited_combined.png"), p_denovo | p_inherited, width = 15, height = 7, dpi = 200)

message("Done. Figures written to: ", out_dir)
