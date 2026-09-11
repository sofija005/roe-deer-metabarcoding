import csv

INPUT_TABLE = "taxonomy/qc/LIB1_ASV_abundance_taxonomy_unfiltered.csv"
INPUT_FASTA = "dada2_out/uncross2/LIB1_ASVs_uncross2.fasta"

OUT_TARGET_TABLE = "taxonomy/processed/LIB1_ASV_abundance_taxonomy_plant_target.csv"
OUT_UNASSIGNED_TABLE = "taxonomy/qc/LIB1_unassigned_ASVs.csv"

OUT_TARGET_FASTA = "taxonomy/processed/LIB1_ASVs_plant_target.fasta"
OUT_UNASSIGNED_FASTA = "taxonomy/qc/LIB1_ASVs_unassigned.fasta"

OUT_SUMMARY = "taxonomy/qc/LIB1_target_taxa_filter_summary.txt"


# ---------------------------------------------------------
# Read taxonomy + abundance table
# ---------------------------------------------------------

target_rows = []
unassigned_rows = []

with open(INPUT_TABLE, newline="", encoding="utf-8-sig") as f:
    reader = csv.DictReader(f)
    fieldnames = reader.fieldnames

    for row in reader:

        phylum = row["Phylum"].strip()

        if phylum == "Tracheophyta":
            target_rows.append(row)

        else:
            unassigned_rows.append(row)


# ---------------------------------------------------------
# Write tables
# ---------------------------------------------------------

with open(OUT_TARGET_TABLE, "w", newline="", encoding="utf-8") as f:
    writer = csv.DictWriter(f, fieldnames=fieldnames)
    writer.writeheader()
    writer.writerows(target_rows)

with open(OUT_UNASSIGNED_TABLE, "w", newline="", encoding="utf-8") as f:
    writer = csv.DictWriter(f, fieldnames=fieldnames)
    writer.writeheader()
    writer.writerows(unassigned_rows)


# ---------------------------------------------------------
# Read FASTA
# ---------------------------------------------------------

sequences = {}

current_id = None
current_seq = []

with open(INPUT_FASTA, encoding="utf-8") as f:
    for line in f:
        line = line.strip()

        if line.startswith(">"):
            if current_id is not None:
                sequences[current_id] = "".join(current_seq)

            current_id = line[1:].split()[0]
            current_seq = []

        else:
            current_seq.append(line)

    if current_id is not None:
        sequences[current_id] = "".join(current_seq)


target_ids = {r["ASV"] for r in target_rows}
unassigned_ids = {r["ASV"] for r in unassigned_rows}


# ---------------------------------------------------------
# Write FASTA
# ---------------------------------------------------------

with open(OUT_TARGET_FASTA, "w") as f:
    for asv in sequences:
        if asv in target_ids:
            f.write(f">{asv}\n{sequences[asv]}\n")

with open(OUT_UNASSIGNED_FASTA, "w") as f:
    for asv in sequences:
        if asv in unassigned_ids:
            f.write(f">{asv}\n{sequences[asv]}\n")


# ---------------------------------------------------------
# Summary
# ---------------------------------------------------------

target_reads = sum(int(r["Total_reads"]) for r in target_rows)
unassigned_reads = sum(int(r["Total_reads"]) for r in unassigned_rows)

total_reads = target_reads + unassigned_reads

low_conf_target = sum(
    1 for r in target_rows
    if r["Assignment_status"] == "Low-confidence"
)

with open(OUT_SUMMARY, "w") as out:

    out.write("=== LIB1 TARGET-TAXA FILTERING ===\n\n")

    out.write(f"Input ASVs: {len(target_rows) + len(unassigned_rows)}\n")
    out.write(f"Input reads: {total_reads}\n\n")

    out.write("Plant target criterion: Phylum = Tracheophyta\n\n")

    out.write(
        f"Plant-target ASVs retained: {len(target_rows)}\n"
    )

    out.write(
        f"Plant-target reads retained: {target_reads} "
        f"({100 * target_reads / total_reads:.2f}%)\n\n"
    )

    out.write(
        f"Unassigned ASVs separated: {len(unassigned_rows)}\n"
    )

    out.write(
        f"Unassigned reads separated: {unassigned_reads} "
        f"({100 * unassigned_reads / total_reads:.2f}%)\n\n"
    )

    out.write(
        f"Low-confidence plant ASVs retained but flagged: "
        f"{low_conf_target}\n"
    )

    out.write(
        "\nNote: unassigned sequences were not interpreted as contaminants. "
        "They were retained separately as unresolved sequences.\n"
    )

print("Target-taxa filtering completed.")
print("Plant target table:", OUT_TARGET_TABLE)
print("Plant target FASTA:", OUT_TARGET_FASTA)
print("Unassigned table:", OUT_UNASSIGNED_TABLE)
print("Summary:", OUT_SUMMARY)
