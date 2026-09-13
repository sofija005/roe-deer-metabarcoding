import csv

table_in = "taxonomy/SINTAX/processed/LIB2_ASV_abundance_taxonomy.csv"
fasta_in = "dada2_out/uncross2/LIB2_ASVs_uncross2.fasta"

table_out = "taxonomy/SINTAX/processed/LIB2_ASV_abundance_taxonomy_plant_target.csv"
fasta_out = "taxonomy/SINTAX/processed/LIB2_ASVs_plant_target.fasta"
summary_out = "taxonomy/SINTAX/qc/LIB2_target_taxa_filter_summary.txt"

target_ids = set()
rows = []

with open(table_in, newline="") as f:
    reader = csv.DictReader(f)
    fields = reader.fieldnames

    for row in reader:
        if row["Phylum"] == "Streptophyta":
            target_ids.add(row["ASV"])
            rows.append(row)

with open(table_out, "w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=fields)
    writer.writeheader()
    writer.writerows(rows)

total_reads = sum(int(r["Total_reads"]) for r in rows)

written = 0
keep = False

with open(fasta_in) as inp, open(fasta_out, "w") as out:
    for line in inp:
        if line.startswith(">"):
            asv = line[1:].strip().split()[0]
            keep = asv in target_ids
            if keep:
                written += 1

        if keep:
            out.write(line)

with open(summary_out, "w") as f:
    f.write("LIB2 TARGET-TAXA FILTER SUMMARY\n")
    f.write("==============================\n\n")
    f.write("Target criterion: Phylum = Streptophyta\n")
    f.write("Reference database: plant-specific Taberlet trnL c/h database\n\n")
    f.write(f"Input ASVs: 313\n")
    f.write(f"Input reads: 236647\n")
    f.write(f"Plant-target ASVs retained: {len(target_ids)}\n")
    f.write(f"Plant-target reads retained: {total_reads}\n")
    f.write(f"Non-target ASVs removed: {313-len(target_ids)}\n\n")
    f.write("Note: Because a plant-specific reference database was used, this step\n")
    f.write("should not be interpreted as an independent test for all possible\n")
    f.write("non-plant contaminants.\n")

print("Plant-target ASVs:", len(target_ids))
print("Plant-target reads:", total_reads)
print("FASTA sequences written:", written)
print("Table:", table_out)
print("FASTA:", fasta_out)
