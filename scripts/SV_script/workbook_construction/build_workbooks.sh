#!/bin/bash
#PBS -l walltime=01:00:00
#PBS -l select=1:ncpus=4:mem=16gb
#PBS -N manta_workbooks
#PBS -o /rds/general/user/rg3825/home/imperial_projects/summer_project/iconik_sv_analysis_public/scripts/SV_script/pbs_logs/
#PBS -e /rds/general/user/rg3825/home/imperial_projects/summer_project/iconik_sv_analysis_public/scripts/SV_script/pbs_logs/

# Stage 3/4: collect the per-family TSVs the qc stage produced and assemble
# the final workbooks - one sheet per family/proband in each.

set -uo pipefail

eval "$(~/miniforge3/bin/conda shell.bash hook)"
conda activate Renv

PROJECT=/rds/general/user/rg3825/home/imperial_projects/summer_project/iconik_sv_analysis_public
SCRIPT_DIR="$PROJECT/scripts/SV_script/workbook_construction"

TRIOS="$PROJECT/data/trios_all37.txt"
OUTDIR="$PROJECT/results/De_novo_calls"
WORKBOOKS_DIR="$PROJECT/results/workbooks"
RSCRIPT=/rds/general/user/rg3825/home/miniforge3/envs/Renv/bin/Rscript

mkdir -p "$WORKBOOKS_DIR"

DENOVO_OUT="$WORKBOOKS_DIR/denovo_R135_all_families.xlsx"
INHERITED_OUT="$WORKBOOKS_DIR/inherited_R135_all_families.xlsx"
ALLGENES_OUT="$WORKBOOKS_DIR/denovo_allgenes_all_families.xlsx"
OVERLAP_OUT="$WORKBOOKS_DIR/denovo_overlap_analysis.xlsx"

DENOVO_TSVS=()
INHERITED_TSVS=()
ALLGENES_TSVS=()

while IFS=$'\t' read -r FAMILY PROBAND_ID PROBAND_SM PROBAND_CRAM MOTHER_ID MOTHER_SM MOTHER_CRAM FATHER_ID FATHER_SM FATHER_CRAM; do

    ANNOTSV_DIR="$OUTDIR/$FAMILY/results/annotsv"
    LOGFILE="$OUTDIR/$FAMILY/${FAMILY}.log"

    if ! grep -q "^Done: $FAMILY\$" "$LOGFILE" 2>/dev/null; then
        echo "WARNING: no confirmed qc completion for $FAMILY (see $LOGFILE) - including any TSVs found anyway"
    fi

    [[ -f "$ANNOTSV_DIR/denovo_R135_green.tsv" ]] && DENOVO_TSVS+=("$ANNOTSV_DIR/denovo_R135_green.tsv")
    [[ -f "$ANNOTSV_DIR/inherited_R135_green.tsv" ]] && INHERITED_TSVS+=("$ANNOTSV_DIR/inherited_R135_green.tsv")
    [[ -f "$ANNOTSV_DIR/denovo_allgenes.tsv" ]] && ALLGENES_TSVS+=("$ANNOTSV_DIR/denovo_allgenes.tsv")

done < "$TRIOS"

if [[ ${#DENOVO_TSVS[@]} -gt 0 ]]; then
    "$RSCRIPT" "$SCRIPT_DIR/build_excel_report.R" --out "$DENOVO_OUT" "${DENOVO_TSVS[@]}"
    echo "De novo workbook written across ${#DENOVO_TSVS[@]} families -> $DENOVO_OUT"
fi

if [[ ${#INHERITED_TSVS[@]} -gt 0 ]]; then
    "$RSCRIPT" "$SCRIPT_DIR/build_excel_report.R" --out "$INHERITED_OUT" "${INHERITED_TSVS[@]}"
    echo "Inherited workbook written across ${#INHERITED_TSVS[@]} families -> $INHERITED_OUT"
fi

if [[ ${#ALLGENES_TSVS[@]} -gt 0 ]]; then
    "$RSCRIPT" "$SCRIPT_DIR/build_excel_report.R" --out "$ALLGENES_OUT" "${ALLGENES_TSVS[@]}"
    echo "Allgenes workbook written across ${#ALLGENES_TSVS[@]} families -> $ALLGENES_OUT"

    "$RSCRIPT" "$SCRIPT_DIR/find_overlapping_denovo.R" --out "$OVERLAP_OUT" "${ALLGENES_TSVS[@]}"
    echo "Overlap workbook written -> $OVERLAP_OUT"
fi
