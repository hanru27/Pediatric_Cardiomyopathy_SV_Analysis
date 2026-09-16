#!/bin/bash
#PBS -l walltime=08:00:00
#PBS -l select=1:ncpus=8:mem=16gb
#PBS -N manta_calling
#PBS -o /rds/general/user/rg3825/home/imperial_projects/summer_project/iconik_sv_analysis_public/scripts/SV_script/pbs_logs/
#PBS -e /rds/general/user/rg3825/home/imperial_projects/summer_project/iconik_sv_analysis_public/scripts/SV_script/pbs_logs/

# Stage 1/4: Manta joint SV calling + de novo scoring, per family.
# Resources are a starting estimate for this stage in isolation - previously
# bundled into one 12h job together with qc, so there's no standalone
# baseline yet; tune after the first run.

set -uo pipefail

source /etc/profile.d/modules.sh
eval "$(~/miniforge3/bin/conda shell.bash hook)"
conda activate manta_env

PROJECT=/rds/general/user/rg3825/home/imperial_projects/summer_project/iconik_sv_analysis_public

# manta_env's Manta build SIGSEGVs during joint calling for some families -
# use the 1.6.0 binary (same as the original paper) for calling, conda env
# only for scoring.
OFFICIAL_MANTA="$PROJECT/tools/manta-1.6.0.centos6_x86_64"
export PATH="$(dirname "$(which python2)"):$PATH"  # official configManta.py needs python2

REF="$PROJECT/reference/Homo_sapiens_assembly38_masked.fasta"
TRIOS="$PROJECT/data/trios_all37.txt"
OUTDIR="$PROJECT/results/De_novo_calls"
THREADS=8

mkdir -p "$OUTDIR"

while IFS=$'\t' read -r FAMILY PROBAND_ID PROBAND_SM PROBAND_CRAM MOTHER_ID MOTHER_SM MOTHER_CRAM FATHER_ID FATHER_SM FATHER_CRAM; do

    echo "Processing family: $FAMILY"
    echo "Proband SM: $PROBAND_SM | Mother SM: $MOTHER_SM | Father SM: $FATHER_SM"

    FAMILY_DIR="$OUTDIR/$FAMILY"
    LOGFILE="$FAMILY_DIR/${FAMILY}.log"
    rm -rf "$FAMILY_DIR"
    mkdir -p "$FAMILY_DIR"

    (
    set -e
    "$OFFICIAL_MANTA/bin/configManta.py" \
        --bam "$PROBAND_CRAM" \
        --bam "$MOTHER_CRAM" \
        --bam "$FATHER_CRAM" \
        --referenceFasta "$REF" \
        --runDir "$FAMILY_DIR"

    "$FAMILY_DIR/runWorkflow.py" -m local -j "$THREADS"

    # denovo_scoring.py uses Python 2 open() and cannot read gzip
    DIPLOID_VCF="$FAMILY_DIR/results/variants/diploidSV.vcf"
    zcat "$FAMILY_DIR/results/variants/diploidSV.vcf.gz" > "$DIPLOID_VCF"

    denovo_scoring.py \
        "$DIPLOID_VCF" \
        "$PROBAND_SM" \
        "$FATHER_SM" \
        "$MOTHER_SM"

    rm -f "$DIPLOID_VCF"

    # left uncompressed - the qc stage reads diploidSV.de_novo.vcf directly;
    # bgzipping here previously made qc silently find nothing (see README)

    echo "Done: $FAMILY"
    ) >> "$LOGFILE" 2>&1 || echo "FAILED: $FAMILY — see $LOGFILE"

done < "$TRIOS"
