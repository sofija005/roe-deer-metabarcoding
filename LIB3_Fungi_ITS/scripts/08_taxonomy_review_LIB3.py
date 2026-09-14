import pandas as pd
import re

TAX = "taxonomy/BOLD_SINTAX/processed/LIB3_SINTAX_taxonomy.csv"
RAW = "taxonomy/BOLD_SINTAX/results/LIB3_SINTAX_all.tsv"
QC = "taxonomy/BOLD_SINTAX/qc"

tax = pd.read_csv(TAX, dtype=str).fillna("")

total = len(tax)

print("=== LIB3 TAXONOMY REVIEW ===")
print("Total ASVs:", total)

# --------------------------------------------------
# 1. Review Family / Genus / Species assignments
# --------------------------------------------------

for rank in ["Family", "Genus", "Species"]:
    assigned = tax[rank].str.strip() != ""
    n = assigned.sum()

    print(
        f"{rank} assigned: {n} "
        f"({100*n/total:.2f}%) | "
        f"unassigned at {rank.lower()}: {total-n}"
    )

    counts = (
        tax.loc[assigned, rank]
        .value_counts()
        .rename_axis(rank)
        .reset_index(name="ASV_count")
    )

    counts.to_csv(
        f"{QC}/LIB3_{rank.lower()}_assignments.csv",
        index=False
    )

# --------------------------------------------------
# 2. Ambiguous taxonomic names
# --------------------------------------------------

ambiguous_pattern = re.compile(
    r"incertae_sedis|uncultured|unidentified|unknown|"
    r"unclassified|environmental|(?:^|_)sp(?:\.|_|$)|"
    r"(?:^|_)cf(?:\.|_|$)|(?:^|_)aff(?:\.|_|$)",
    re.IGNORECASE
)

def is_ambiguous(row):
    text = " ".join([
        row["Family"],
        row["Genus"],
        row["Species"]
    ])
    return bool(ambiguous_pattern.search(text))

ambiguous = tax[tax.apply(is_ambiguous, axis=1)].copy()

ambiguous.to_csv(
    f"{QC}/LIB3_ambiguous_assignments.csv",
    index=False
)

# --------------------------------------------------
# 3. Low-confidence proposed assignments
#    Look directly at raw SINTAX output
# --------------------------------------------------

raw = pd.read_csv(
    RAW,
    sep="\t",
    header=None,
    names=["ASV_ID", "Raw_assignment", "Strand", "Accepted_taxonomy"],
    dtype=str
).fillna("")

rank_names = {
    "k": "Kingdom",
    "p": "Phylum",
    "c": "Class",
    "o": "Order",
    "f": "Family",
    "g": "Genus",
    "s": "Species"
}

low_rows = []

pattern = re.compile(
    r"([kpcofgs]):([^,()]+)\(([0-9.]+)\)"
)

for _, row in raw.iterrows():
    for prefix, taxon, conf in pattern.findall(row["Raw_assignment"]):

        conf = float(conf)

        if conf < 0.80:
            low_rows.append({
                "ASV_ID": row["ASV_ID"],
                "Rank": rank_names[prefix],
                "Proposed_taxon": taxon,
                "Confidence": conf
            })

low = pd.DataFrame(low_rows)

if len(low):
    low.to_csv(
        f"{QC}/LIB3_low_confidence_assignments.csv",
        index=False
    )

    low_asvs = low["ASV_ID"].nunique()
else:
    low_asvs = 0

# --------------------------------------------------
# 4. Summary
# --------------------------------------------------

print("\nAmbiguous ASVs:", len(ambiguous))
print("ASVs with at least one proposed rank below 0.80:", low_asvs)

print("\n=== TOP FAMILIES ===")
print(
    tax.loc[tax["Family"] != "", "Family"]
    .value_counts()
    .head(15)
    .to_string()
)

print("\n=== TOP GENERA ===")
print(
    tax.loc[tax["Genus"] != "", "Genus"]
    .value_counts()
    .head(15)
    .to_string()
)

print("\n=== TOP SPECIES ===")
print(
    tax.loc[tax["Species"] != "", "Species"]
    .value_counts()
    .head(15)
    .to_string()
)

print("\nQC files written to:", QC)
