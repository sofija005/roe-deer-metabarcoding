import csv

INFILE = "taxonomy/qc/LIB1_ASV_abundance_taxonomy_unfiltered.csv"

LOW = "taxonomy/qc/LIB1_BOLD_low_confidence.csv"
AMB_GENUS = "taxonomy/qc/LIB1_BOLD_ambiguous_genus.csv"
AMB_SPECIES = "taxonomy/qc/LIB1_BOLD_ambiguous_species.csv"
SPECIES = "taxonomy/qc/LIB1_BOLD_species_assignments.csv"
SUMMARY = "taxonomy/qc/LIB1_BOLD_final_QC_summary.txt"

with open(INFILE, newline="", encoding="utf-8-sig") as f:
    reader = csv.DictReader(f)
    fields = reader.fieldnames
    rows = list(reader)

low = [r for r in rows if r["Assignment_status"] == "Low-confidence"]
amb_g = [r for r in rows if r["Ambiguous_genus"] == "YES"]
amb_s = [r for r in rows if r["Ambiguous_species"] == "YES"]
species = [r for r in rows if r["Tax_Rank"] == "SPECIES"]

plant = [r for r in rows if r["Phylum"] == "Tracheophyta"]
unassigned = [r for r in rows if not r["Phylum"].strip()]

def total_reads(data):
    return sum(int(r["Total_reads"]) for r in data)

def write_csv(path, data):
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(data)

write_csv(LOW, low)
write_csv(AMB_GENUS, amb_g)
write_csv(AMB_SPECIES, amb_s)
write_csv(SPECIES, species)

with open(SUMMARY, "w", encoding="utf-8") as out:
    out.write("=== LIB1 BOLD / BLAST FINAL QC ===\n\n")

    out.write(f"Total ASVs: {len(rows)}\n")
    out.write(f"Total reads: {total_reads(rows)}\n\n")

    out.write(f"Plant-target ASVs: {len(plant)}\n")
    out.write(f"Plant-target reads: {total_reads(plant)}\n\n")

    out.write(f"Unassigned ASVs: {len(unassigned)}\n")
    out.write(f"Unassigned reads: {total_reads(unassigned)}\n\n")

    out.write(f"Low-confidence ASVs: {len(low)}\n")
    out.write(f"Low-confidence reads: {total_reads(low)}\n\n")

    out.write(f"Ambiguous genus ASVs: {len(amb_g)}\n")
    out.write(f"Ambiguous genus reads: {total_reads(amb_g)}\n\n")

    out.write(f"Ambiguous species ASVs: {len(amb_s)}\n")
    out.write(f"Ambiguous species reads: {total_reads(amb_s)}\n\n")

    out.write(f"Species-level assignments: {len(species)}\n")
    out.write(f"Species-level reads: {total_reads(species)}\n\n")

    out.write(
        "Unassigned sequences were retained separately as unresolved "
        "and were not interpreted as contaminants.\n"
    )

print("Final BOLD QC completed.")
print("Summary:", SUMMARY)
