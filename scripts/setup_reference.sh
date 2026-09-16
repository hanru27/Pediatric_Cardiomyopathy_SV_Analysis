#!/bin/bash
# Downloads the GRCh38 DRAGEN masked reference (public, no auth needed) and
# builds the .fai/.dict/BWA index this pipeline expects in reference/.
# Source: gs://gcp-public-data--broad-references/hg38/v0/dragen_reference/
# (Broad's public reference bucket - see GATK resource bundle docs).
#
# Requires: wget, samtools, bwa on PATH.
#
# Usage:
#   ./setup_reference.sh [reference_dir]   # defaults to ../reference

set -euo pipefail

REF_DIR=${1:-"$(dirname "$0")/../reference"}
BUCKET_URL=https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/dragen_reference
FASTA=Homo_sapiens_assembly38_masked.fasta

mkdir -p "$REF_DIR"
cd "$REF_DIR"

if [[ ! -f "$FASTA" ]]; then
    echo "Downloading $FASTA ..."
    wget -q --show-progress "$BUCKET_URL/$FASTA" -O "$FASTA"
fi

echo "Building .fai ..."
samtools faidx "$FASTA"

echo "Building .dict ..."
samtools dict "$FASTA" > "${FASTA%.fasta}.dict"

echo "Building BWA index (this takes a while on a 3GB genome) ..."
bwa index "$FASTA"

echo "Done. Reference ready in $REF_DIR"
