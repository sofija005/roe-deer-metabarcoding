# LIB2 – Plants trnL

This directory contains the bioinformatic processing workflow for the
plant trnL metabarcoding library of roe deer (*Capreolus capreolus*).

## Workflow

1. Primer removal with Cutadapt
2. Quality filtering with DADA2
3. Error learning and denoising
4. Paired-end merging and ASV table construction
5. Chimera removal
6. UNCROSS2 tag-jump filtering
7. Final sequence and ASV tracking

## Results

The final LIB2 dataset contained:

- 301,578 input reads
- 255,980 quality-filtered reads
- 247,640 merged reads
- 236,650 non-chimeric reads
- 236,647 reads after UNCROSS2 filtering
- 368 raw ASVs
- 313 final ASVs

Contaminant assessment was not performed because suitable negative controls
or associated metadata were not available.

## Directories

- `scripts/` – scripts used for each processing step
- `results/` – summary tables, final ASV abundance table, and final ASV FASTA file
