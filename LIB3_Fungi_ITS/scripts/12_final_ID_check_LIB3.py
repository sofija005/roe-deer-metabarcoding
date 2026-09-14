import pandas as pd

FASTA = "taxonomy/BOLD_SINTAX/filtered/LIB3_ASVs_target_Fungi.fasta"
ABUND = "taxonomy/BOLD_SINTAX/filtered/LIB3_ASV_table_target_Fungi.csv"
TAX = "taxonomy/BOLD_SINTAX/processed/LIB3_SINTAX_taxonomy.csv"
FINAL = "taxonomy/BOLD_SINTAX/final/LIB3_final_ASV_abundance_taxonomy.csv"

OUT = "taxonomy/BOLD_SINTAX/qc/LIB3_final_ID_check.txt"

# ----------------------------
# FASTA IDs
# ----------------------------
fasta_ids = []

with open(FASTA) as f:
    for line in f:
        if line.startswith(">"):
            fasta_ids.append(line[1:].strip().split()[0])

# ----------------------------
# Table IDs
# ----------------------------
def read_ids(path):
    df = pd.read_csv(path, dtype=str)

    if "ASV_ID" not in df.columns:
        df = df.rename(columns={df.columns[0]: "ASV_ID"})

    return df["ASV_ID"].tolist()

abund_ids = read_ids(ABUND)
tax_ids = read_ids(TAX)
final_ids = read_ids(FINAL)

# Convert to sets for comparison
fasta_set = set(fasta_ids)
abund_set = set(abund_ids)
tax_set = set(tax_ids)
final_set = set(final_ids)

# ----------------------------
# Write summary
# ----------------------------
with open(OUT, "w") as f:
    f.write("=== LIB3 FINAL ASV ID CONSISTENCY CHECK ===\n\n")

    f.write(f"FASTA IDs: {len(fasta_ids)}\n")
    f.write(f"Abundance table IDs: {len(abund_ids)}\n")
    f.write(f"Taxonomy table IDs: {len(tax_ids)}\n")
    f.write(f"Final merged table IDs: {len(final_ids)}\n\n")

    f.write(f"FASTA == abundance: {fasta_set == abund_set}\n")
    f.write(f"FASTA == taxonomy: {fasta_set == tax_set}\n")
    f.write(f"FASTA == final table: {fasta_set == final_set}\n")
    f.write(f"Abundance == taxonomy: {abund_set == tax_set}\n")
    f.write(f"Abundance == final table: {abund_set == final_set}\n\n")

    f.write(f"IDs missing from FASTA: {len(final_set - fasta_set)}\n")
    f.write(f"IDs missing from abundance table: {len(final_set - abund_set)}\n")
    f.write(f"IDs missing from taxonomy table: {len(final_set - tax_set)}\n")
    f.write(f"IDs missing from final table: {len(fasta_set - final_set)}\n\n")

    f.write(f"Duplicate FASTA IDs: {len(fasta_ids) - len(fasta_set)}\n")
    f.write(f"Duplicate abundance IDs: {len(abund_ids) - len(abund_set)}\n")
    f.write(f"Duplicate taxonomy IDs: {len(tax_ids) - len(tax_set)}\n")
    f.write(f"Duplicate final-table IDs: {len(final_ids) - len(final_set)}\n")

print(open(OUT).read())
