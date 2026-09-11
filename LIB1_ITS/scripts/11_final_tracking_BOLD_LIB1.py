import csv
from collections import defaultdict

TRACKING = "/home/sofija/INTERNSHIP/roe-deer-metabarcoding/LIB1_ITS/results/tracking/LIB1_final_sequence_tracking.csv"

TAX = "taxonomy/qc/LIB1_ASV_abundance_taxonomy_unfiltered.csv"

OUT = "taxonomy/qc/LIB1_final_sequence_tracking_BOLD.csv"

plant = defaultdict(int)
unassigned = defaultdict(int)

# Read ASV abundance + taxonomy table
with open(TAX, newline="", encoding="utf-8-sig") as f:
    reader = csv.DictReader(f)

    samples = [x for x in reader.fieldnames if x.startswith("LME")]

    for row in reader:

        is_plant = row["Phylum"].strip() == "Tracheophyta"

        for sample in samples:

            n = int(row[sample])

            if is_plant:
                plant[sample] += n
            else:
                unassigned[sample] += n


# Read original sequence tracking
with open(TRACKING, newline="", encoding="utf-8-sig") as f:
    reader = csv.DictReader(f)

    rows = list(reader)
    old_fields = reader.fieldnames


new_fields = old_fields + [
    "BOLD_plant_target",
    "BOLD_unassigned"
]


for row in rows:

    sample = row["Sample"]

    row["BOLD_plant_target"] = plant[sample]
    row["BOLD_unassigned"] = unassigned[sample]

    uncross = int(row["UNCROSS2_filtered"])

    # Everything after UNCROSS2 must end up
    # either as plant-target or unassigned
    if plant[sample] + unassigned[sample] != uncross:

        raise ValueError(
            f"Tracking mismatch for {sample}: "
            f"{plant[sample]} + {unassigned[sample]} != {uncross}"
        )


with open(OUT, "w", newline="", encoding="utf-8") as f:

    writer = csv.DictWriter(f, fieldnames=new_fields)

    writer.writeheader()
    writer.writerows(rows)


print("Tracking PASS")
print("Plant target:", sum(plant.values()))
print("Unassigned:", sum(unassigned.values()))
print("Total:", sum(plant.values()) + sum(unassigned.values()))
