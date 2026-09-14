import pandas as pd
from pathlib import Path

ABUND = Path("taxonomy/BOLD_SINTAX/filtered/LIB3_ASV_table_target_Fungi.csv")
TAX = Path("taxonomy/BOLD_SINTAX/processed/LIB3_SINTAX_taxonomy.csv")

OUTDIR = Path("taxonomy/BOLD_SINTAX/final")
OUTDIR.mkdir(parents=True, exist_ok=True)

OUT = OUTDIR / "LIB3_final_ASV_abundance_taxonomy.csv"

abund = pd.read_csv(ABUND)
tax = pd.read_csv(TAX, dtype=str).fillna("")

# Detect/standardize ASV ID column
if "ASV_ID" not in abund.columns:
    abund = abund.rename(columns={abund.columns[0]: "ASV_ID"})

# Keep taxonomy fields we want in final table
tax_cols = [
    "ASV_ID",
    "Kingdom",
    "Kingdom_confidence",
    "Phylum",
    "Phylum_confidence",
    "Class",
    "Class_confidence",
    "Order",
    "Order_confidence",
    "Family",
    "Family_confidence",
    "Genus",
    "Genus_confidence",
    "Species",
    "Species_confidence",
    "Strand",
    "Deepest_rank",
    "Deepest_taxon",
    "Accepted_taxonomy"
]

tax = tax[tax_cols]

merged = abund.merge(
    tax,
    on="ASV_ID",
    how="left",
    validate="one_to_one"
)

missing_tax = merged["Kingdom"].isna().sum()

merged.to_csv(OUT, index=False)

print("=== FINAL LIB3 ABUNDANCE + TAXONOMY TABLE ===")
print("ASVs in abundance table:", len(abund))
print("ASVs in taxonomy table:", len(tax))
print("ASVs in final merged table:", len(merged))
print("ASVs missing taxonomy:", missing_tax)
print("Output:", OUT)
