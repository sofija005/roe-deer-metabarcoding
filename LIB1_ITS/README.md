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
8. BOLD taxonomic assignment using the PLANT:PUBLIC reference library.
9. Merge BOLD classification and combined-hit outputs across all batches.
10. Link taxonomy assignments to the UNCROSS2-filtered ASV abundance table.
11. Perform taxonomy QC:
    - review species, genus and family assignments,
    - flag unassigned, ambiguous and low-confidence assignments,
    - inspect unexpected taxa,
    - evaluate sample-level assignment quality.
12. Filter to the biological plant target (Tracheophyta).
13. Generate the final plant-target ASV abundance + taxonomy table.
14. Generate the final plant-target FASTA and separate unresolved ASVs.
15. Perform final sequence tracking and FASTA/table ASV ID consistency checks.

## Directories

- `scripts/` – scripts used for each processing step
- `results/` – summary tables and final tracking results
