# LIB2 – Plants trnL

This directory contains the bioinformatic processing workflow for the plant trnL
metabarcoding library of roe deer (*Capreolus capreolus*).

## Workflow

1. Primer removal with Cutadapt
2. Quality filtering with DADA2
3. Error learning and denoising
4. Paired-end merging and ASV table construction
5. Chimera removal
6. UNCROSS2 tag-jump filtering
7. Final sequence and ASV tracking
7. Taxonomic assignment using SINTAX with a marker-specific Taberlet trnL c/h reference database.
8. Apply a SINTAX confidence threshold of 0.80.
9. Link accepted taxonomy assignments to the UNCROSS2-filtered ASV abundance table.
10. Review species-, genus-, family- and higher-rank assignments.
11. Flag unresolved lower-rank assignments and inspect possible non-target taxa.
12. Retain plant-target sequences classified within Streptophyta.
13. Generate the final plant-target ASV abundance + taxonomy table and FASTA.
14. Perform final taxonomy-aware sequence tracking.
15. Verify consistency of ASV identifiers between the final FASTA and abundance table.

## Directories

- `scripts/` – scripts used for each processing step
- `results/denoising/` – DADA2 denoising outputs
- `results/merged/` – merged paired-end output
- `results/asv_tables/` – raw, chimera-filtered and UNCROSS2-filtered ASV tables
- `results/fasta/` – final ASV sequences in FASTA format
- `results/tracking/` – sequence/ASV tracking and filtering summaries
