#!/bin/bash
#PBS -l walltime=04:00:00
#PBS -l select=1:ncpus=16:mem=32gb
#PBS -N manta_samplots
#PBS -o /rds/general/user/rg3825/home/imperial_projects/summer_project/iconik_sv_analysis_public/scripts/SV_script/pbs_logs/
#PBS -e /rds/general/user/rg3825/home/imperial_projects/summer_project/iconik_sv_analysis_public/scripts/SV_script/pbs_logs/

# Stage 4/4: samplot trio images for every R135 green-gene candidate the qc stage produced - one PNG per unique SV, plus a self-contained index.html gallery per family/label.

set -uo pipefail

PROJECT=/rds/general/user/rg3825/home/imperial_projects/summer_project/iconik_sv_analysis_public
SCRIPT_DIR="$PROJECT/scripts/figure_scripts/samplots"

TRIOS="$PROJECT/data/trios_all37.txt"
REF="$PROJECT/reference/Homo_sapiens_assembly38_masked.fasta"
OUTDIR="$PROJECT/results/De_novo_calls"
VIS_OUTDIR="$PROJECT/results/variant_visualisations"

# samplot needs a newer python than manta_env/Renv so make sure to call by full path instead of activating over manta_env
SAMPLOT_PYTHON=/rds/general/user/rg3825/home/miniforge3/envs/samplot_env/bin/python3
SAMPLOT_BIN=/rds/general/user/rg3825/home/miniforge3/envs/samplot_env/bin/samplot

MAX_PARALLEL=16

while IFS=$'\t' read -r FAMILY PROBAND_ID PROBAND_SM PROBAND_CRAM MOTHER_ID MOTHER_SM MOTHER_CRAM FATHER_ID FATHER_SM FATHER_CRAM; do

    FAMILY_DIR="$OUTDIR/$FAMILY"
    LOGFILE="$FAMILY_DIR/${FAMILY}.log"
    ANNOTSV_DIR="$FAMILY_DIR/results/annotsv"
    DENOVO_TSV="$ANNOTSV_DIR/denovo_R135_green.tsv"
    INHERITED_TSV="$ANNOTSV_DIR/inherited_R135_green.tsv"
    DENOVO_VIS_DIR="$VIS_OUTDIR/$FAMILY/denovo_visualisations"
    # kept misspelled to match the folder already on disk
    INHERITED_VIS_DIR="$VIS_OUTDIR/$FAMILY/inherited_visulaisations"

    if [[ ! -f "$DENOVO_TSV" && ! -f "$INHERITED_TSV" ]]; then
        echo "WARNING: no R135 green TSVs for $FAMILY, skipping samplot"
        continue
    fi

    mkdir -p "$DENOVO_VIS_DIR" "$INHERITED_VIS_DIR"

    (
    set -e
    # de novo samplot
    "$SAMPLOT_PYTHON" "$SCRIPT_DIR/samplot_trio.py" \
        --tsv "$DENOVO_TSV" --label denovo --family "$FAMILY" \
        --proband-cram "$PROBAND_CRAM" --proband-id "$PROBAND_SM" \
        --mother-cram  "$MOTHER_CRAM"  --mother-id  "$MOTHER_SM" \
        --father-cram  "$FATHER_CRAM"  --father-id  "$FATHER_SM" \
        --reference "$REF" \
        --outdir "$DENOVO_VIS_DIR" \
        --samplot-bin "$SAMPLOT_BIN"

    # inherited samplot
    "$SAMPLOT_PYTHON" "$SCRIPT_DIR/samplot_trio.py" \
        --tsv "$INHERITED_TSV" --label inherited --family "$FAMILY" \
        --proband-cram "$PROBAND_CRAM" --proband-id "$PROBAND_SM" \
        --mother-cram  "$MOTHER_CRAM"  --mother-id  "$MOTHER_SM" \
        --father-cram  "$FATHER_CRAM"  --father-id  "$FATHER_SM" \
        --reference "$REF" \
        --outdir "$INHERITED_VIS_DIR" \
        --samplot-bin "$SAMPLOT_BIN"

    # de novo gallery
    "$SAMPLOT_PYTHON" "$SCRIPT_DIR/build_samplot_gallery.py" \
        --pngdir "$DENOVO_VIS_DIR" --family "$FAMILY" --label denovo
    # inherited gallery
    "$SAMPLOT_PYTHON" "$SCRIPT_DIR/build_samplot_gallery.py" \
        --pngdir "$INHERITED_VIS_DIR" --family "$FAMILY" --label inherited

    echo "Done: samplot $FAMILY"
    ) >> "$LOGFILE" 2>&1 &

    while (( $(jobs -rp | wc -l) >= MAX_PARALLEL )); do
        wait -n
    done

done < "$TRIOS"

wait
