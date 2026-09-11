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

## Directories

- `scripts/` – scripts used for each processing step
- `results/denoising/` – DADA2 denoising outputs
- `results/merged/` – merged paired-end output
- `results/asv_tables/` – raw, chimera-filtered and UNCROSS2-filtered ASV tables
- `results/fasta/` – final ASV sequences in FASTA format
- `results/tracking/` – sequence/ASV tracking and filtering summaries
