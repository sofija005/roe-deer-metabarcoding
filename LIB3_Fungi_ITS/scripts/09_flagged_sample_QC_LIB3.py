import pandas as pd
from pathlib import Path

BASE = Path("taxonomy/BOLD_SINTAX")
QC = BASE / "qc"

abund_path = "dada2_out/uncross2/LIB3_ASV_table_uncross2.csv"
tax_path = BASE / "processed/LIB3_SINTAX_taxonomy.csv"

amb_path = QC / "LIB3_ambiguous_assignments.csv"
low_path = QC / "LIB3_low_confidence_assignments.csv"

abund = pd.read_csv(abund_path)
tax = pd.read_csv(tax_path, dtype=str).fillna("")

# Detect ASV ID column
if "ASV_ID" in abund.columns:
    id_col = "ASV_ID"
else:
    id_col = abund.columns[0]
    abund = abund.rename(columns={id_col: "ASV_ID"})
    id_col = "ASV_ID"

sample_cols = [c for c in abund.columns if c != id_col]

# Make sure abundance columns are numeric
for c in sample_cols:
    abund[c] = pd.to_numeric(abund[c], errors="coerce").fillna(0)

def summarize_flagged(flag_file, label):
    if not Path(flag_file).exists():
        print(f"{label}: file not found")
        return

    flags = pd.read_csv(flag_file, dtype=str).fillna("")

    ids = set(flags["ASV_ID"])
    x = abund[abund["ASV_ID"].isin(ids)].copy()

    if x.empty:
        print(f"{label}: no matching ASVs")
        return

    rows = []

    for _, row in x.iterrows():
        asv = row["ASV_ID"]

        present = {
            s: int(row[s])
            for s in sample_cols
            if row[s] > 0
        }

        rows.append({
            "ASV_ID": asv,
            "Total_reads": sum(present.values()),
            "Samples_present": len(present),
            "Sample_list": ";".join(present.keys()),
            "Sample_read_counts": ";".join(
                f"{s}:{n}" for s, n in present.items()
            )
        })

    out = pd.DataFrame(rows)

    out.to_csv(
        QC / f"LIB3_{label}_sample_association.csv",
        index=False
    )

    print(f"\n=== {label.upper()} SAMPLE ASSOCIATION ===")
    print("Flagged ASVs:", len(out))

    sample_summary = {}

    for sample in sample_cols:
        n = (x[sample] > 0).sum()
        reads = int(x[sample].sum())

        if n > 0:
            sample_summary[sample] = (n, reads)

    summary_df = pd.DataFrame(
        [
            {
                "Sample": s,
                "Flagged_ASVs": n,
                "Flagged_reads": r
            }
            for s, (n, r) in sample_summary.items()
        ]
    ).sort_values(
        ["Flagged_ASVs", "Flagged_reads"],
        ascending=False
    )

    summary_df.to_csv(
        QC / f"LIB3_{label}_by_sample.csv",
        index=False
    )

    print(summary_df.head(20).to_string(index=False))


summarize_flagged(amb_path, "ambiguous")
summarize_flagged(low_path, "low_confidence")

print("\nQC files written to:", QC)
