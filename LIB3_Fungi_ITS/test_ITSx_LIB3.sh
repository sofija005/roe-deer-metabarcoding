#!/bin/bash

set -euo pipefail

BASE="/home/sofija/INTERNSHIP/samples/Capreolus_capreolus/Fungi_ITS"

INPUT="$BASE/test_dada2_out/test_ASVs.fasta"
OUTDIR="$BASE/test_itsx_out"

mkdir -p "$OUTDIR"

ITSx \
    -i "$INPUT" \
    -o "$OUTDIR/fungi_test" \
    -t F

echo "ITSx test completed successfully."
