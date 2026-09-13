import csv

tracking_in = "/home/sofija/INTERNSHIP/roe-deer-metabarcoding/LIB2_trnL/results/tracking/LIB2_final_sequence_tracking.csv"
taxonomy_table = "taxonomy/SINTAX/processed/LIB2_ASV_abundance_taxonomy_plant_target.csv"

tracking_out = "taxonomy/SINTAX/qc/LIB2_final_sequence_tracking_SINTAX.csv"

with open(taxonomy_table, newline="") as f:
    reader = csv.DictReader(f)
    tax_rows = list(reader)
    sample_cols = [
        c for c in reader.fieldnames
        if c.startswith("LME")
    ]

target_reads = {
    sample: sum(int(r[sample]) for r in tax_rows)
    for sample in sample_cols
}

with open(tracking_in, newline="") as f:
    reader = csv.DictReader(f)
    tracking = list(reader)
    old_fields = reader.fieldnames

new_fields = old_fields + [
    "Plant_target_after_taxonomy",
    "Non_target_removed_after_taxonomy"
]

all_pass = True

for row in tracking:
    sample = row["Sample"]

    uncross = int(row["UNCROSS2_filtered"])
    plant = target_reads.get(sample, 0)
    removed = uncross - plant

    row["Plant_target_after_taxonomy"] = plant
    row["Non_target_removed_after_taxonomy"] = removed

    if plant + removed != uncross:
        all_pass = False

with open(tracking_out, "w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=new_fields)
    writer.writeheader()
    writer.writerows(tracking)

print("Tracking file:", tracking_out)
print("Samples:", len(tracking))
print("Plant-target reads:", sum(target_reads.values()))
print("RESULT:", "PASS" if all_pass else "FAIL")
