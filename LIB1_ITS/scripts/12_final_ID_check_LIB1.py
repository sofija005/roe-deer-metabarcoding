import csv

TARGET_TABLE = "taxonomy/processed/LIB1_ASV_abundance_taxonomy_plant_target.csv"
TARGET_FASTA = "taxonomy/processed/LIB1_ASVs_plant_target.fasta"

UNASSIGNED_TABLE = "taxonomy/qc/LIB1_unassigned_ASVs.csv"
UNASSIGNED_FASTA = "taxonomy/qc/LIB1_ASVs_unassigned.fasta"

OUT = "taxonomy/qc/LIB1_final_FASTA_table_ID_check.txt"

def table_ids(path):
    with open(path, newline="", encoding="utf-8-sig") as f:
        return [r["ASV"] for r in csv.DictReader(f)]

def fasta_ids(path):
    ids = []
    with open(path, encoding="utf-8") as f:
        for line in f:
            if line.startswith(">"):
                ids.append(line[1:].strip().split()[0])
    return ids

target_t = table_ids(TARGET_TABLE)
target_f = fasta_ids(TARGET_FASTA)

unass_t = table_ids(UNASSIGNED_TABLE)
unass_f = fasta_ids(UNASSIGNED_FASTA)

target_missing_fasta = set(target_t) - set(target_f)
target_missing_table = set(target_f) - set(target_t)

unass_missing_fasta = set(unass_t) - set(unass_f)
unass_missing_table = set(unass_f) - set(unass_t)

with open(OUT, "w") as out:
    out.write("=== LIB1 FINAL FASTA / TABLE ID CHECK ===\n\n")

    out.write("PLANT TARGET\n")
    out.write(f"Table ASVs: {len(target_t)}\n")
    out.write(f"FASTA ASVs: {len(target_f)}\n")
    out.write(f"Missing from FASTA: {len(target_missing_fasta)}\n")
    out.write(f"Missing from table: {len(target_missing_table)}\n\n")

    out.write("UNASSIGNED\n")
    out.write(f"Table ASVs: {len(unass_t)}\n")
    out.write(f"FASTA ASVs: {len(unass_f)}\n")
    out.write(f"Missing from FASTA: {len(unass_missing_fasta)}\n")
    out.write(f"Missing from table: {len(unass_missing_table)}\n\n")

    if (
        not target_missing_fasta
        and not target_missing_table
        and not unass_missing_fasta
        and not unass_missing_table
    ):
        out.write("RESULT: PASS\n")
    else:
        out.write("RESULT: FAIL\n")

print("Done.")
