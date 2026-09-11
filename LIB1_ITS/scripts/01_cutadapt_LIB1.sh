#!/bin/bash

# =========================
# Plants ITS - LIB1
# Primer removal with Cutadapt
# =========================

FWD_PRIMER="ATGCGATACTTGGTGTGAAT"
REV_PRIMER="GACGCTTCTCCAGACTACAAT"

INPUT_DIR="merged_raw"
OUTPUT_DIR="primersCut_out"
LOG_DIR="logs"

mkdir -p "$OUTPUT_DIR"
mkdir -p "$LOG_DIR"

for R1 in "$INPUT_DIR"/*_R1.fq.gz
do
    sample=$(basename "$R1" "_R1.fq.gz")
    R2="$INPUT_DIR/${sample}_R2.fq.gz"

    echo "Processing ${sample}..."

    cutadapt \
        -g "^${FWD_PRIMER}" \
        -G "^${REV_PRIMER}" \
        -o "$OUTPUT_DIR/${sample}_R1.fq.gz" \
        -p "$OUTPUT_DIR/${sample}_R2.fq.gz" \
        --minimum-length 32 \
        --pair-filter=both \
        --cores=0 \
        "$R1" "$R2" \
        2>&1 | tee "$LOG_DIR/${sample}_cutadapt.log"

done

echo "Cutadapt finished for LIB1."
