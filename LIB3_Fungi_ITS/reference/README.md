# BOLD fungal ITS reference database

The large BOLD reference files are not stored in this repository.

Reference source:
- BOLD Systems public fungal records
- Kingdom: Fungi
- Markers retained: ITS, ITS1 and ITS2

The BOLD TSV export is converted into a VSEARCH SINTAX-compatible FASTA database using:

`../scripts/07a_prepare_BOLD_reference_LIB3.py`

Expected local database filename:

`BOLD_Fungi_ITS_SINTAX.fasta`

Taxonomic classification is performed with VSEARCH SINTAX using a confidence cutoff of 0.80 and both strands.
