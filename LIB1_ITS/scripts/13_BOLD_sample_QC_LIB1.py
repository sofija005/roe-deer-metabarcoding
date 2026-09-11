import csv

INFILE = "taxonomy/qc/LIB1_ASV_abundance_taxonomy_unfiltered.csv"

OUT = "taxonomy/qc/LIB1_BOLD_sample_QC.csv"


with open(INFILE, newline="", encoding="utf-8-sig") as f:

    reader = csv.DictReader(f)

    samples = [
        x for x in reader.fieldnames
        if x.startswith("LME")
    ]

    rows = list(reader)


results = []


for sample in samples:

    total = 0
    plant = 0
    unassigned = 0
    low_conf = 0
    ambiguous_genus = 0
    ambiguous_species = 0

    for row in rows:

        n = int(row[sample])

        total += n

        if row["Phylum"].strip() == "Tracheophyta":
            plant += n
        else:
            unassigned += n

        if row["Assignment_status"] == "Low-confidence":
            low_conf += n

        if row["Ambiguous_genus"] == "YES":
            ambiguous_genus += n

        if row["Ambiguous_species"] == "YES":
            ambiguous_species += n


    results.append({

        "Sample": sample,

        "Total_UNCROSS2_reads": total,

        "Plant_target_reads": plant,

        "Plant_target_percent":
            round(100 * plant / total, 2) if total else 0,

        "Unassigned_reads": unassigned,

        "Unassigned_percent":
            round(100 * unassigned / total, 2) if total else 0,

        "Low_confidence_reads": low_conf,

        "Ambiguous_genus_reads": ambiguous_genus,

        "Ambiguous_species_reads": ambiguous_species
    })


with open(OUT, "w", newline="", encoding="utf-8") as f:

    writer = csv.DictWriter(
        f,
        fieldnames=results[0].keys()
    )

    writer.writeheader()
    writer.writerows(results)


print("=== SAMPLE QC ===")

for r in results:

    print(
        r["Sample"],
        "plant =", str(r["Plant_target_percent"]) + "%",
        "unassigned =", str(r["Unassigned_percent"]) + "%"
    )

print("\nSaved:", OUT)
