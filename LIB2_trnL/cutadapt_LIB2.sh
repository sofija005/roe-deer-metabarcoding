#!/bin/bash

# =========================
# Plants trnLc-trnLh - LIB2
# Primer removal with Cutadapt
# =========================

FWD_PRIMER="CGAAATCGGTAGACGCTACG"
REV_PRIMER="CCATTGAGTCTCTGCACCTATC"

OUTPUT_DIR="primersCut_out"
LOG_DIR="logs"

mkdir -p "$OUTPUT_DIR"
mkdir -p "$LOG_DIR"

for d in LME*
do
    sample=$(basename "$d")

    R1=$(find "$d" -maxdepth 1 -type f -name "*_1.fq.gz")
    R2=$(find "$d" -maxdepth 1 -type f -name "*_2.fq.gz")

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

echo "Cutadapt finished for LIB2."
