import csv
from collections import Counter, defaultdict

asv_table = "dada2_out/uncross2/LIB2_ASV_table_uncross2.csv"
tax_file = "taxonomy/SINTAX/processed/LIB2_SINTAX_taxonomy_QC.csv"

out_file = "taxonomy/SINTAX/processed/LIB2_ASV_abundance_taxonomy.csv"
summary_file = "taxonomy/SINTAX/qc/LIB2_taxonomy_abundance_QC_summary.txt"

# -----------------------------
# Read taxonomy
# -----------------------------
taxonomy = {}

with open(tax_file, newline="") as f:
    reader = csv.DictReader(f)
    for row in reader:
        taxonomy[row["ASV"]] = row

# -----------------------------
# Read abundance table
# -----------------------------
with open(asv_table, newline="") as f:
    reader = csv.DictReader(f)
    sample_columns = [x for x in reader.fieldnames if x != "ASV"]
    abundance_rows = list(reader)

abundance_ids = {r["ASV"] for r in abundance_rows}
taxonomy_ids = set(taxonomy)

missing_taxonomy = sorted(abundance_ids - taxonomy_ids)
extra_taxonomy = sorted(taxonomy_ids - abundance_ids)

rank_asvs = Counter()
rank_reads = Counter()

phylum_asvs = Counter()
phylum_reads = Counter()

family_asvs = Counter()
family_reads = Counter()

genus_asvs = Counter()
genus_reads = Counter()

species_asvs = Counter()
species_reads = Counter()

output_rows = []

total_reads = 0

for row in abundance_rows:
    asv = row["ASV"]
    total = sum(int(row[s]) for s in sample_columns)
    total_reads += total

    tax = taxonomy.get(asv, {})

    accepted_rank = tax.get("Accepted_rank", "Unassigned")
    phylum = tax.get("Phylum", "") or "Unassigned"
    family = tax.get("Family", "") or "Unassigned"
    genus = tax.get("Genus", "") or "Unassigned"
    species = tax.get("Species", "") or "Unassigned"

    rank_asvs[accepted_rank] += 1
    rank_reads[accepted_rank] += total

    phylum_asvs[phylum] += 1
    phylum_reads[phylum] += total

    family_asvs[family] += 1
    family_reads[family] += total

    genus_asvs[genus] += 1
    genus_reads[genus] += total

    species_asvs[species] += 1
    species_reads[species] += total

    out = {
        "ASV": asv,
        "Total_reads": total,
    }

    for s in sample_columns:
        out[s] = row[s]

    for col in [
        "Kingdom", "Kingdom_confidence",
        "Phylum", "Phylum_confidence",
        "Class", "Class_confidence",
        "Order", "Order_confidence",
        "Family", "Family_confidence",
        "Genus", "Genus_confidence",
        "Species", "Species_confidence",
        "Accepted_rank", "Accepted_confidence",
        "Strand"
    ]:
        out[col] = tax.get(col, "")

    output_rows.append(out)

fields = (
    ["ASV", "Total_reads"]
    + sample_columns
    + [
        "Kingdom", "Kingdom_confidence",
        "Phylum", "Phylum_confidence",
        "Class", "Class_confidence",
        "Order", "Order_confidence",
        "Family", "Family_confidence",
        "Genus", "Genus_confidence",
        "Species", "Species_confidence",
        "Accepted_rank", "Accepted_confidence",
        "Strand"
    ]
)

with open(out_file, "w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=fields)
    writer.writeheader()
    writer.writerows(output_rows)

with open(summary_file, "w") as f:
    f.write("LIB2 TAXONOMY + ABUNDANCE QC SUMMARY\n")
    f.write("====================================\n\n")

    f.write(f"Total ASVs in abundance table: {len(abundance_rows)}\n")
    f.write(f"Total reads: {total_reads}\n")
    f.write(f"Missing taxonomy IDs: {len(missing_taxonomy)}\n")
    f.write(f"Extra taxonomy IDs: {len(extra_taxonomy)}\n\n")

    f.write("LOWEST ACCEPTED TAXONOMIC RANK\n")
    for rank, n in rank_asvs.most_common():
        reads = rank_reads[rank]
        pct = reads / total_reads * 100 if total_reads else 0
        f.write(f"{rank}: {n} ASVs, {reads} reads ({pct:.2f}%)\n")

    f.write("\nPHYLA\n")
    for taxon, n in phylum_asvs.most_common():
        reads = phylum_reads[taxon]
        pct = reads / total_reads * 100 if total_reads else 0
        f.write(f"{taxon}: {n} ASVs, {reads} reads ({pct:.2f}%)\n")

    f.write("\nTOP FAMILIES BY READ ABUNDANCE\n")
    for taxon, reads in family_reads.most_common(15):
        f.write(f"{taxon}: {family_asvs[taxon]} ASVs, {reads} reads\n")

    f.write("\nTOP GENERA BY READ ABUNDANCE\n")
    for taxon, reads in genus_reads.most_common(15):
        f.write(f"{taxon}: {genus_asvs[taxon]} ASVs, {reads} reads\n")

    f.write("\nSPECIES-LEVEL ASSIGNMENTS\n")
    for taxon, reads in species_reads.most_common():
        if taxon != "Unassigned":
            f.write(f"{taxon}: {species_asvs[taxon]} ASVs, {reads} reads\n")

print("Output table:", out_file)
print("Summary:", summary_file)
print("ASVs:", len(abundance_rows))
print("Reads:", total_reads)
print("Missing taxonomy IDs:", len(missing_taxonomy))
print("Extra taxonomy IDs:", len(extra_taxonomy))
