import pandas as pd
from pathlib import Path

TAX = Path("taxonomy/BOLD_SINTAX/processed/LIB3_SINTAX_taxonomy.csv")
ABUND = Path("dada2_out/uncross2/LIB3_ASV_table_uncross2.csv")
FASTA = Path("dada2_out/uncross2/LIB3_ASVs_uncross2.fasta")

OUT = Path("taxonomy/BOLD_SINTAX/filtered")
QC = Path("taxonomy/BOLD_SINTAX/qc")

OUT.mkdir(parents=True, exist_ok=True)

# --------------------------------------------------
# 1. Read taxonomy
# --------------------------------------------------

tax = pd.read_csv(TAX, dtype=str).fillna("")

target_ids = set(
    tax.loc[tax["Kingdom"].str.strip() == "Fungi", "ASV_ID"]
)

unassigned = tax[tax["Kingdom"].str.strip() == ""].copy()

non_target = tax[
    (tax["Kingdom"].str.strip() != "") &
    (tax["Kingdom"].str.strip() != "Fungi")
].copy()

# --------------------------------------------------
# 2. Read abundance table
# --------------------------------------------------

abund = pd.read_csv(ABUND)

# Detect orientation safely
if "ASV_ID" in abund.columns:
    id_col = "ASV_ID"

elif abund.iloc[:, 0].astype(str).str.startswith("ASV").any():
    abund = abund.rename(columns={abund.columns[0]: "ASV_ID"})
    id_col = "ASV_ID"

else:
    # ASVs may be columns instead of rows
    asv_columns = [c for c in abund.columns if str(c).startswith("ASV")]

    if asv_columns:
        abund = abund.set_index(abund.columns[0]).T.reset_index()
        abund = abund.rename(columns={"index": "ASV_ID"})
        id_col = "ASV_ID"
    else:
        raise ValueError("Could not detect ASV orientation in abundance table.")

sample_cols = [c for c in abund.columns if c != id_col]

for c in sample_cols:
    abund[c] = pd.to_numeric(abund[c], errors="coerce").fillna(0)

asvs_before = len(abund)
reads_before = int(abund[sample_cols].sum().sum())

filtered = abund[abund[id_col].astype(str).isin(target_ids)].copy()

asvs_after = len(filtered)
reads_after = int(filtered[sample_cols].sum().sum())

filtered.to_csv(
    OUT / "LIB3_ASV_table_target_Fungi.csv",
    index=False
)

# --------------------------------------------------
# 3. Save non-target / unassigned lists
# --------------------------------------------------

non_target.to_csv(
    QC / "LIB3_non_target_ASVs.csv",
    index=False
)

unassigned.to_csv(
    QC / "LIB3_unassigned_ASVs.csv",
    index=False
)

# --------------------------------------------------
# 4. Filter FASTA
# --------------------------------------------------

kept_fasta = 0

with open(FASTA) as fin, \
     open(OUT / "LIB3_ASVs_target_Fungi.fasta", "w") as fout:

    keep = False

    for line in fin:
        if line.startswith(">"):
            seq_id = line[1:].strip().split()[0]
            keep = seq_id in target_ids

            if keep:
                kept_fasta += 1
                fout.write(line)

        elif keep:
            fout.write(line)

# --------------------------------------------------
# 5. Summary
# --------------------------------------------------

summary = QC / "LIB3_target_taxa_filter_summary.txt"

with open(summary, "w") as f:
    f.write("=== LIB3 TARGET-TAXA FILTERING ===\n")
    f.write("Target kingdom: Fungi\n\n")

    f.write(f"ASVs before filtering: {asvs_before}\n")
    f.write(f"ASVs after filtering: {asvs_after}\n")
    f.write(
        f"ASV retention: {100*asvs_after/asvs_before:.2f}%\n\n"
    )

    f.write(f"Reads before filtering: {reads_before}\n")
    f.write(f"Reads after filtering: {reads_after}\n")
    f.write(
        f"Read retention: {100*reads_after/reads_before:.2f}%\n\n"
    )

    f.write(f"Non-target ASVs: {len(non_target)}\n")
    f.write(f"Unassigned ASVs: {len(unassigned)}\n")
    f.write(f"FASTA sequences retained: {kept_fasta}\n")

print("=== LIB3 TARGET-TAXA FILTERING COMPLETED ===")
print("ASVs before:", asvs_before)
print("ASVs after:", asvs_after)
print("Reads before:", reads_before)
print("Reads after:", reads_after)
print("Non-target ASVs:", len(non_target))
print("Unassigned ASVs:", len(unassigned))
print("FASTA sequences retained:", kept_fasta)
