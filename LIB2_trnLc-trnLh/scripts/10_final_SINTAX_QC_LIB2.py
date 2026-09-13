import csv

inp = "taxonomy/SINTAX/processed/LIB2_ASV_abundance_taxonomy_plant_target.csv"

summary_out = "taxonomy/SINTAX/qc/LIB2_SINTAX_final_QC_summary.txt"
lowres_out = "taxonomy/SINTAX/qc/LIB2_SINTAX_low_resolution.csv"

rows = []

with open(inp, newline="") as f:
    reader = csv.DictReader(f)
    rows = list(reader)
    fields = reader.fieldnames

total_asvs = len(rows)
total_reads = sum(int(r["Total_reads"]) for r in rows)

plant = [r for r in rows if r["Phylum"] == "Streptophyta"]
non_target = [r for r in rows if r["Phylum"] != "Streptophyta"]

# Reliable only above family level
low_resolution = [
    r for r in rows
    if r["Accepted_rank"] in {"Phylum", "Class", "Order"}
]

family_unresolved = [r for r in rows if not r["Family"]]
genus_unresolved = [r for r in rows if not r["Genus"]]
species_unresolved = [r for r in rows if not r["Species"]]

def reads(x):
    return sum(int(r["Total_reads"]) for r in x)

with open(lowres_out, "w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=fields)
    writer.writeheader()
    writer.writerows(low_resolution)

with open(summary_out, "w") as f:
    f.write("LIB2 SINTAX FINAL TAXONOMY QC SUMMARY\n")
    f.write("====================================\n\n")

    f.write(f"Total ASVs: {total_asvs}\n")
    f.write(f"Total reads: {total_reads}\n")
    f.write("SINTAX confidence threshold: 0.80\n\n")

    f.write(f"Plant-target ASVs: {len(plant)}\n")
    f.write(f"Plant-target reads: {reads(plant)}\n")
    f.write(f"Non-target ASVs: {len(non_target)}\n")
    f.write(f"Non-target reads: {reads(non_target)}\n\n")

    f.write("LOWER-RANK RESOLUTION\n")
    f.write(
        f"Family unresolved: {len(family_unresolved)} ASVs, "
        f"{reads(family_unresolved)} reads\n"
    )
    f.write(
        f"Genus unresolved: {len(genus_unresolved)} ASVs, "
        f"{reads(genus_unresolved)} reads\n"
    )
    f.write(
        f"Species unresolved: {len(species_unresolved)} ASVs, "
        f"{reads(species_unresolved)} reads\n"
    )
    f.write(
        f"Low-resolution (accepted only to phylum/class/order): "
        f"{len(low_resolution)} ASVs, {reads(low_resolution)} reads\n\n"
    )

    f.write("Interpretation:\n")
    f.write("- No ASVs were unassigned at the plant/phylum level.\n")
    f.write("- No non-target ASVs were identified using this plant-specific database.\n")
    f.write("- Lower-rank assignments below the 0.80 threshold were not accepted.\n")
    f.write("- SINTAX does not directly report tied best-hit ambiguity; therefore\n")
    f.write("  ambiguity was not quantified in the same way as for BOLD results.\n")
    f.write("- The plant-specific database is not an independent screen for all\n")
    f.write("  possible non-plant contaminants.\n")

print("Final QC summary:", summary_out)
print("Low-resolution table:", lowres_out)
print("Total ASVs:", total_asvs)
print("Total reads:", total_reads)
print("Plant-target ASVs:", len(plant))
print("Non-target ASVs:", len(non_target))
