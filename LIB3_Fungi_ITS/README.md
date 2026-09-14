# LIB3 – Fungi ITS

This directory contains the bioinformatic processing workflow for the fungal ITS
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

## Taxonomic assignment and final QC

Taxonomic assignment was performed on the UNCROSS2-filtered LIB3 Fungi_ITS ASVs.

### Reference database

Public fungal records were downloaded from BOLD Systems. Records belonging to the ITS, ITS1 and ITS2 markers were retained and converted into a VSEARCH SINTAX-compatible FASTA reference database.

The large BOLD reference files are not stored in this repository. Instructions for preparing the local reference database are provided in `reference/README.md`.

### Taxonomy workflow

1. `07a_prepare_BOLD_reference_LIB3.py`
   - prepares the local BOLD fungal ITS reference database for SINTAX.

2. `07b_run_SINTAX_LIB3.sh`
   - performs taxonomic classification with VSEARCH SINTAX.
   - confidence threshold: 0.80
   - both strands are considered.

3. `07_SINTAX_taxonomy_QC_LIB3.py`
   - converts the SINTAX output into a structured taxonomy table.
   - records the deepest accepted taxonomic rank.

4. `08_taxonomy_review_LIB3.py`
   - reviews family-, genus- and species-level assignments.
   - flags ambiguous and low-confidence assignments.

5. `09_flagged_sample_QC_LIB3.py`
   - evaluates the distribution of ambiguous and low-confidence ASVs across samples.

6. `10_target_taxa_filter_LIB3.py`
   - retains ASVs consistent with the fungal biological target.

7. `11_final_abundance_taxonomy_LIB3.py`
   - combines the final ASV abundance table with the accepted taxonomy.

8. `12_final_ID_check_LIB3.py`
   - verifies consistency of ASV IDs between FASTA, abundance and taxonomy files.

### Final outputs

Final outputs include:

- structured SINTAX taxonomy table
- target-Fungi ASV abundance table
- target-Fungi FASTA
- final ASV abundance + taxonomy table
- ambiguous and low-confidence assignment summaries
- sample-level taxonomy QC
- final sequence-tracking summary
- final ASV ID consistency check
