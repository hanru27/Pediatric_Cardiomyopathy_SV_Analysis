#!/bin/bash
#PBS -l walltime=06:00:00
#PBS -l select=1:ncpus=16:mem=32gb
#PBS -N manta_qc
#PBS -o /rds/general/user/rg3825/home/imperial_projects/summer_project/iconik_sv_analysis/scripts/SV_script/pbs_logs/
#PBS -e /rds/general/user/rg3825/home/imperial_projects/summer_project/iconik_sv_analysis/scripts/SV_script/pbs_logs/

# Stage 2/4: tier-1 filtering + de novo/inherited split, then AnnotSV +
# SpliceAI/GREEN-DB against the R135 panel and genome-wide. Produces
# per-family TSVs; workbook_construction/ assembles the final workbooks.
# SpliceAI on inherited candidates is the dominant cost (~30-50min/family).

set -uo pipefail

source /etc/profile.d/modules.sh
eval "$(~/miniforge3/bin/conda shell.bash hook)"
conda activate manta_env
module load BCFtools/1.22-GCC-14.2.0
module load BEDTools/2.31.1-GCC-14.3.0

PROJECT=/rds/general/user/rg3825/home/imperial_projects/summer_project/iconik_sv_analysis
SCRIPT_DIR="$PROJECT/scripts/SV_script/qc"

REF="$PROJECT/reference/Homo_sapiens_assembly38_masked.fasta"
TRIOS="$PROJECT/data/trios_all37.txt"
PEDIGREE="$PROJECT/data/iconik_trioPedigree.csv"
OUTDIR="$PROJECT/results/De_novo_calls"

RSCRIPT=/rds/general/user/rg3825/home/miniforge3/envs/Renv/bin/Rscript
ANNOTSV_BIN="$PROJECT/tools/AnnotSV-3.5.10/bin/AnnotSV"
ANNOTSV_ANNOTATIONS="$PROJECT/tools/AnnotSV-3.5.10/bin/AnnotSV_annotations"
GREEN_GENES="$PROJECT/data/R135/R135_green_genes.tsv"

# SpliceAI/GREEN-DB run in their own conda env (TensorFlow, python 3.9) -
# called by full interpreter path so manta_env stays active for the rest of this script
SPLICEAI_PYTHON=/rds/general/user/rg3825/home/miniforge3/envs/spliceai_env/bin/python3
ANNOTATE_SPLICE_GREENDB="$SCRIPT_DIR/annotate_splice_greendb.py"
GREENVARAN_BIN="$PROJECT/tools/GREEN-VARAN/bin/greenvaran"
GREENDB_BED="$PROJECT/tools/GREEN-VARAN/resources/GRCh38/GRCh38_GREEN-DB.bed.gz"
GREENDB_SCHEMA="$PROJECT/tools/GREEN-VARAN/config/greendb_schema_v2.5.json"

MAX_PARALLEL=16

####===========
# step 1: tier-1 filtering (PASS, no BND, 50bp-25kb, proband carries it), then de novo/inherited split

while IFS=$'\t' read -r FAMILY PROBAND_ID PROBAND_SM PROBAND_CRAM MOTHER_ID MOTHER_SM MOTHER_CRAM FATHER_ID FATHER_SM FATHER_CRAM; do

    FAMILY_DIR="$OUTDIR/$FAMILY"
    LOGFILE="$FAMILY_DIR/${FAMILY}.log"
    DENOVO_VCF="$FAMILY_DIR/results/variants/diploidSV.de_novo.vcf"

    (
    set -e
    TIER1_INPUT="$DENOVO_VCF"
    TIER1_OUTPUT="$FAMILY_DIR/results/variants/diploidSV.tier1.vcf.gz"

    bcftools view \
        -f PASS \
        -e 'SVTYPE="BND" || ABS(SVLEN) < 50 || ABS(SVLEN) > 25000' \
        "$TIER1_INPUT" \
        -O z -o "$TIER1_OUTPUT"

    bcftools index "$TIER1_OUTPUT"

    TIER1_COUNT=$(bcftools view -H "$TIER1_OUTPUT" | wc -l)
    echo "Tier-1 filtered variants: $TIER1_COUNT"

    # Summary table
    TIER1_SUMMARY="$FAMILY_DIR/results/variants/diploidSV.tier1_summary.tsv"
    echo -e "CHROM\tPOS\tSVTYPE\tSVLEN\tFILTER" > "$TIER1_SUMMARY"
    bcftools query \
        -f '%CHROM\t%POS\t%INFO/SVTYPE\t%INFO/SVLEN\t%FILTER\n' \
        "$TIER1_OUTPUT" >> "$TIER1_SUMMARY"

    # proband carries it (GT=alt, FT=PASS) - not gated on DQ=60 (Manta's conflict-only flag),
    # that hides ordinary inherited SVs; de novo vs inherited is decided from raw GT below.
    # -s reorders samples so GT[0]/FT[0]=proband, GT[1]=mother, GT[2]=father throughout
    ALL_CANDIDATES="$FAMILY_DIR/results/variants/diploidSV.tier1_all_candidates.vcf.gz"

    bcftools view \
        -s "$PROBAND_SM,$MOTHER_SM,$FATHER_SM" \
        -i 'GT[0]="alt" && FMT/FT[0]="PASS"' \
        "$TIER1_OUTPUT" \
        -O z -o "$ALL_CANDIDATES"

    bcftools index -f "$ALL_CANDIDATES"

    ALL_COUNT=$(bcftools view -H "$ALL_CANDIDATES" | wc -l)
    echo "Tier-1 candidates where proband carries the variant: $ALL_COUNT"

    # de novo = neither parent carries it; inherited = at least one does. Exception: non-PAR chrX,
    # male proband only checks mum (dad can't contribute an X to a son).
    # trios.txt IDs sometimes have a trailing _D<n> the pedigree file doesn't - strip before matching
    PROBAND_ID_NO_SUFFIX=$(echo "$PROBAND_ID" | sed -E 's/_D[0-9]+$//')
    PROBAND_SEX=$(awk -F',' -v id="$PROBAND_ID_NO_SUFFIX" '$2==id {print $3}' "$PEDIGREE")

    if [[ "$PROBAND_SEX" != "Male" && "$PROBAND_SEX" != "Female" ]]; then
        echo "WARNING: proband sex not found in pedigree for $FAMILY ($PROBAND_ID) - chrX hemizygous rule NOT applied (dad's genotype will count on chrX too), check $PEDIGREE manually"
    fi

    DENOVO_CANDIDATES="$FAMILY_DIR/results/variants/diploidSV.tier1_denovo_candidates.vcf.gz"
    INHERITED_CANDIDATES="$FAMILY_DIR/results/variants/diploidSV.tier1_inherited_candidates.vcf.gz"

    if [[ "$PROBAND_SEX" == "Male" ]]; then
        # GRCh38 PAR1/PAR2 spelled out in full - bcftools can't parse a negated parenthesised sub-expression
        INHERITED_EXPR='(CHROM="chrX" && (POS<10001 || POS>2781479) && (POS<155701383 || POS>156030895) && GT[1]!="ref") || ((CHROM!="chrX" || (CHROM="chrX" && POS>=10001 && POS<=2781479) || (CHROM="chrX" && POS>=155701383 && POS<=156030895)) && (GT[1]!="ref" || GT[2]!="ref"))'
    else
        # female proband (or unknown) - chrX is normal diploid, either parent counts
        INHERITED_EXPR='GT[1]!="ref" || GT[2]!="ref"'
    fi

    bcftools view -i "$INHERITED_EXPR" "$ALL_CANDIDATES" -O z -o "$INHERITED_CANDIDATES"
    bcftools index -f "$INHERITED_CANDIDATES"

    bcftools view -e "$INHERITED_EXPR" "$ALL_CANDIDATES" -O z -o "$DENOVO_CANDIDATES"
    bcftools index -f "$DENOVO_CANDIDATES"

    DENOVO_COUNT=$(bcftools view -H "$DENOVO_CANDIDATES" | wc -l)
    INHERITED_COUNT=$(bcftools view -H "$INHERITED_CANDIDATES" | wc -l)
    echo "Tier-1 de novo candidates: $DENOVO_COUNT"
    echo "Tier-1 inherited candidates: $INHERITED_COUNT"

    echo "Done: $FAMILY"
    ) >> "$LOGFILE" 2>&1 || echo "FAILED: $FAMILY — see $LOGFILE"

done < "$TRIOS"

####===========
# step 2: AnnotSV annotation + R135 green gene filtering (de novo + inherited)

# runs AnnotSV + R135 filtering on one candidate VCF - called once for de novo, once for inherited
annotate_and_filter() {
    local VCF="$1"
    local LABEL="$2"
    local GT_TSV="$ANNOTSV_DIR/${LABEL}_genotypes.tsv"
    FILTERED_TSV="$ANNOTSV_DIR/${LABEL}_R135_green.tsv"

    if [[ ! -f "$VCF" ]]; then
        echo "WARNING: $VCF not found, skipping $LABEL for $FAMILY"
        return
    fi

    # column names carry role + sample ID, so genotype ownership is obvious without cross-referencing trios.txt
    printf 'CHROM\tPOS\tSVTYPE\tProband_%s_GT\tProband_%s_GQ\tProband_%s_FT\tProband_%s_PR\tProband_%s_SR\tMother_%s_GT\tMother_%s_GQ\tMother_%s_FT\tMother_%s_PR\tMother_%s_SR\tFather_%s_GT\tFather_%s_GQ\tFather_%s_FT\tFather_%s_PR\tFather_%s_SR\n' \
        "$PROBAND_SM" "$PROBAND_SM" "$PROBAND_SM" "$PROBAND_SM" "$PROBAND_SM" \
        "$MOTHER_SM"  "$MOTHER_SM"  "$MOTHER_SM"  "$MOTHER_SM"  "$MOTHER_SM" \
        "$FATHER_SM"  "$FATHER_SM"  "$FATHER_SM"  "$FATHER_SM"  "$FATHER_SM" \
        > "$GT_TSV"
    bcftools query \
        -s "${PROBAND_SM},${MOTHER_SM},${FATHER_SM}" \
        -f '%CHROM\t%POS\t%INFO/SVTYPE\t[%GT\t%GQ\t%FT\t%PR\t%SR\t]\n' \
        "$VCF" | sed 's/\t$//' >> "$GT_TSV"

    echo "$LABEL genotypes extracted: $(tail -n +2 "$GT_TSV" | wc -l) variants"

    "$ANNOTSV_BIN" \
        -SVinputFile    "$VCF" \
        -outputFile     "$ANNOTSV_DIR/${LABEL}_annotated" \
        -genomeBuild    GRCh38 \
        -SVminSize      50 \
        -annotationsDir "$ANNOTSV_ANNOTATIONS"

    # SpliceAI (breakpoint proxy scoring) + GREEN-DB, one row per SV, joined below by Manta SV ID
    SPLICE_GREENDB_TSV="$ANNOTSV_DIR/${LABEL}_splice_greendb.tsv"
    CUDA_VISIBLE_DEVICES="" "$SPLICEAI_PYTHON" "$ANNOTATE_SPLICE_GREENDB" \
        --vcf              "$VCF" \
        --ref              "$REF" \
        --greenvaran-bin   "$GREENVARAN_BIN" \
        --greendb          "$GREENDB_BED" \
        --greendb-schema   "$GREENDB_SCHEMA" \
        --out              "$SPLICE_GREENDB_TSV"

    # filters to R135 green genes, attaches genotypes + SpliceAI/GREEN-DB (see filter_annotsv_panel.R)
    "$RSCRIPT" "$SCRIPT_DIR/filter_annotsv_panel.R" \
        --annotsv         "$ANNOTSV_DIR/${LABEL}_annotated.tsv" \
        --genes           "$GREEN_GENES" \
        --gt              "$GT_TSV" \
        --family          "$FAMILY" \
        --proband         "$PROBAND_SM" \
        --pedigree        "$PEDIGREE" \
        --proband_id      "$PROBAND_ID" \
        --mother_id       "$MOTHER_ID" \
        --father_id       "$FATHER_ID" \
        --splice_greendb  "$SPLICE_GREENDB_TSV" \
        --out             "$FILTERED_TSV"
}

while IFS=$'\t' read -r FAMILY PROBAND_ID PROBAND_SM PROBAND_CRAM MOTHER_ID MOTHER_SM MOTHER_CRAM FATHER_ID FATHER_SM FATHER_CRAM; do

    FAMILY_DIR="$OUTDIR/$FAMILY"
    LOGFILE="$FAMILY_DIR/${FAMILY}.log"
    DENOVO_CANDIDATES="$FAMILY_DIR/results/variants/diploidSV.tier1_denovo_candidates.vcf.gz"
    INHERITED_CANDIDATES="$FAMILY_DIR/results/variants/diploidSV.tier1_inherited_candidates.vcf.gz"
    ANNOTSV_DIR="$FAMILY_DIR/results/annotsv"

    if [[ ! -f "$DENOVO_CANDIDATES" ]]; then
        echo "WARNING: $DENOVO_CANDIDATES not found, skipping $FAMILY"
        continue
    fi

    mkdir -p "$ANNOTSV_DIR"

    # run this family in the background so the next one can start immediately
    (
    set -e
    annotate_and_filter "$DENOVO_CANDIDATES" "denovo"
    annotate_and_filter "$INHERITED_CANDIDATES" "inherited"
    echo "Done: $FAMILY"
    ) >> "$LOGFILE" 2>&1 &

    # throttle to MAX_PARALLEL concurrent families so we don't oversubscribe the node
    while (( $(jobs -rp | wc -l) >= MAX_PARALLEL )); do
        wait -n
    done

done < "$TRIOS"

wait

####===========
# step 3: genome-wide (no panel) de novo filtering, reusing step 2's AnnotSV/genotype TSVs

while IFS=$'\t' read -r FAMILY PROBAND_ID PROBAND_SM PROBAND_CRAM MOTHER_ID MOTHER_SM MOTHER_CRAM FATHER_ID FATHER_SM FATHER_CRAM; do

    FAMILY_DIR="$OUTDIR/$FAMILY"
    LOGFILE="$FAMILY_DIR/${FAMILY}.log"
    ANNOTSV_DIR="$FAMILY_DIR/results/annotsv"
    ANNOTSV_TSV="$ANNOTSV_DIR/denovo_annotated.tsv"
    GT_TSV="$ANNOTSV_DIR/denovo_genotypes.tsv"
    ALLGENES_TSV="$ANNOTSV_DIR/denovo_allgenes.tsv"

    if [[ ! -f "$ANNOTSV_TSV" ]]; then
        echo "WARNING: $ANNOTSV_TSV not found, skipping allgenes for $FAMILY"
        continue
    fi

    (
    set -e
    "$RSCRIPT" "$SCRIPT_DIR/filter_annotsv_allgenes.R" \
        --annotsv    "$ANNOTSV_TSV" \
        --gt         "$GT_TSV" \
        --family     "$FAMILY" \
        --proband    "$PROBAND_SM" \
        --pedigree   "$PEDIGREE" \
        --proband_id "$PROBAND_ID" \
        --mother_id  "$MOTHER_ID" \
        --father_id  "$FATHER_ID" \
        --out        "$ALLGENES_TSV"
    echo "Done: allgenes $FAMILY"
    ) >> "$LOGFILE" 2>&1 &

    while (( $(jobs -rp | wc -l) >= MAX_PARALLEL )); do
        wait -n
    done

done < "$TRIOS"

wait
