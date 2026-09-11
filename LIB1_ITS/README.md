# LIB1 – Plants ITS

This directory contains the bioinformatic processing workflow for the
plant ITS metabarcoding library of roe deer (*Capreolus capreolus*).

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
- `results/` – summary tables and final tracking results
