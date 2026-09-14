#!/bin/bash

set -euo pipefail

QUERY="results/fasta/LIB3_ASVs_uncross2.fasta"
DB="reference/BOLD_Fungi_ITS_SINTAX.fasta"
OUT="results/taxonomy/LIB3_SINTAX_all.tsv"

mkdir -p results/taxonomy

vsearch \
  --sintax "$QUERY" \
  --db "$DB" \
  --tabbedout "$OUT" \
  --sintax_cutoff 0.80 \
  --strand both \
  --threads 8
