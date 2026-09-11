import csv
from collections import defaultdict

INFILE = "taxonomy/qc/LIB1_ASV_abundance_taxonomy_unfiltered.csv"

OUT = "taxonomy/qc/LIB1_taxa_for_manual_review.csv"


summary = defaultdict(lambda: {
    "ASVs": 0,
    "Reads": 0
})


with open(INFILE, newline="", encoding="utf-8-sig") as f:

    reader = csv.DictReader(f)

    for row in reader:

        if row["Species"].strip():

            rank = "SPECIES"
            taxon = row["Species"].strip()

        elif row["Genus"].strip():

            rank = "GENUS"
            taxon = row["Genus"].strip()

        elif row["Family"].strip():

            rank = "FAMILY"
            taxon = row["Family"].strip()

        elif row["Order"].strip():

            rank = "ORDER"
            taxon = row["Order"].strip()

        elif row["Class"].strip():

            rank = "CLASS"
            taxon = row["Class"].strip()

        else:

            rank = "UNASSIGNED"
            taxon = "Unassigned"


        key = (rank, taxon)

        summary[key]["ASVs"] += 1
        summary[key]["Reads"] += int(row["Total_reads"])


rows = []

for (rank, taxon), values in summary.items():

    rows.append({

        "Taxonomic_rank": rank,
        "Taxon": taxon,
        "ASVs": values["ASVs"],
        "Reads": values["Reads"],
        "Review_status": "",
        "Notes": ""
    })


rows.sort(
    key=lambda x: x["Reads"],
    reverse=True
)


with open(OUT, "w", newline="", encoding="utf-8") as f:

    writer = csv.DictWriter(
        f,
        fieldnames=rows[0].keys()
    )

    writer.writeheader()
    writer.writerows(rows)


print("Taxa review table created:", OUT)

print("\nTop assignments:")

for r in rows[:25]:

    print(
        r["Taxonomic_rank"],
        r["Taxon"],
        "-",
        r["ASVs"],
        "ASVs,",
        r["Reads"],
        "reads"
    )
