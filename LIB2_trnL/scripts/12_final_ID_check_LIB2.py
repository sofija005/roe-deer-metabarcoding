import csv

table = "taxonomy/SINTAX/processed/LIB2_ASV_abundance_taxonomy_plant_target.csv"
fasta = "taxonomy/SINTAX/processed/LIB2_ASVs_plant_target.fasta"
out = "taxonomy/SINTAX/qc/LIB2_final_FASTA_table_ID_check.txt"

with open(table, newline="") as f:
    table_ids = {r["ASV"] for r in csv.DictReader(f)}

fasta_ids = set()

with open(fasta) as f:
    for line in f:
        if line.startswith(">"):
            fasta_ids.add(line[1:].strip().split()[0])

missing_fasta = table_ids - fasta_ids
missing_table = fasta_ids - table_ids

result = (
    "PASS"
    if not missing_fasta and not missing_table
    else "FAIL"
)

with open(out, "w") as f:
    f.write("LIB2 FINAL FASTA/TABLE ID CHECK\n")
    f.write("===============================\n\n")
    f.write(f"Table ASVs: {len(table_ids)}\n")
    f.write(f"FASTA ASVs: {len(fasta_ids)}\n")
    f.write(f"Missing from FASTA: {len(missing_fasta)}\n")
    f.write(f"Missing from table: {len(missing_table)}\n")
    f.write(f"RESULT: {result}\n")

print(open(out).read())
