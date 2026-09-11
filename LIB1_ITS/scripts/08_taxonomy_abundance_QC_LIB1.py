import csv
from collections import Counter, defaultdict

ABUND_FILE = "dada2_out/uncross2/LIB1_ASV_table_uncross2.csv"
TAX_FILE   = "taxonomy/qc/LIB1_BOLD_taxonomy_QC.csv"

OUT_TABLE = "taxonomy/qc/LIB1_ASV_abundance_taxonomy_unfiltered.csv"
OUT_SUMMARY = "taxonomy/qc/LIB1_taxonomy_abundance_QC_summary.txt"

# ---------------------------------------------------------
# Taxonomy
# ---------------------------------------------------------

taxonomy = {}

with open(TAX_FILE, newline="", encoding="utf-8-sig") as f:
    reader = csv.DictReader(f)
    tax_fields = reader.fieldnames

    for row in reader:
        taxonomy[row["ASV"]] = row

# ---------------------------------------------------------
# Abundance table
# ---------------------------------------------------------

results = []

with open(ABUND_FILE, newline="", encoding="utf-8-sig") as f:
    reader = csv.DictReader(f)

    sample_names = [x for x in reader.fieldnames if x != "ASV"]

    for row in reader:

        asv = row["ASV"]

        counts = {}
        total_reads = 0

        for sample in sample_names:
            value = int(float(row[sample]))
            counts[sample] = value
            total_reads += value

        tax = taxonomy.get(asv)

        if tax is None:
            print("WARNING: no taxonomy for", asv)
            continue

        combined = {
            "ASV": asv,
            "Total_reads": total_reads
        }

        for sample in sample_names:
            combined[sample] = counts[sample]

        for field in tax_fields:
            if field != "ASV":
                combined[field] = tax[field]

        results.append(combined)

# ---------------------------------------------------------
# ID consistency
# ---------------------------------------------------------

abundance_ids = {r["ASV"] for r in results}
taxonomy_ids = set(taxonomy)

missing_taxonomy = abundance_ids - taxonomy_ids
extra_taxonomy = taxonomy_ids - abundance_ids

# ---------------------------------------------------------
# Write joined table
# ---------------------------------------------------------

fieldnames = list(results[0].keys())

with open(OUT_TABLE, "w", newline="", encoding="utf-8") as f:
    writer = csv.DictWriter(f, fieldnames=fieldnames)
    writer.writeheader()
    writer.writerows(results)

# ---------------------------------------------------------
# Summaries
# ---------------------------------------------------------

def summarize(field):

    asv_counts = Counter()
    read_counts = Counter()

    for r in results:

        value = r[field].strip() if r[field] else "Unassigned"

        asv_counts[value] += 1
        read_counts[value] += int(r["Total_reads"])

    return asv_counts, read_counts


status_asv, status_reads = summarize("Assignment_status")
phylum_asv, phylum_reads = summarize("Phylum")
family_asv, family_reads = summarize("Family")
genus_asv, genus_reads = summarize("Genus")
species_asv, species_reads = summarize("Species")

total_reads = sum(int(r["Total_reads"]) for r in results)

ambiguous_genus_asvs = [
    r for r in results
    if r["Ambiguous_genus"] == "YES"
]

ambiguous_species_asvs = [
    r for r in results
    if r["Ambiguous_species"] == "YES"
]

with open(OUT_SUMMARY, "w", encoding="utf-8") as out:

    out.write("=== LIB1 TAXONOMY + ABUNDANCE QC ===\n\n")

    out.write(f"Total ASVs: {len(results)}\n")
    out.write(f"Total reads: {total_reads}\n")
    out.write(f"Missing taxonomy IDs: {len(missing_taxonomy)}\n")
    out.write(f"Extra taxonomy IDs: {len(extra_taxonomy)}\n")

    out.write("\n=== ASSIGNMENT STATUS ===\n")

    for status in status_asv:
        reads = status_reads[status]

        out.write(
            f"{status}: "
            f"{status_asv[status]} ASVs, "
            f"{reads} reads "
            f"({100 * reads / total_reads:.2f}%)\n"
        )

    out.write("\n=== PHYLUM ===\n")

    for name, n in phylum_asv.most_common():
        reads = phylum_reads[name]

        out.write(
            f"{name}: {n} ASVs, "
            f"{reads} reads "
            f"({100 * reads / total_reads:.2f}%)\n"
        )

    out.write("\n=== TOP 20 FAMILIES BY READS ===\n")

    for name, reads in family_reads.most_common(20):
        out.write(
            f"{name}: {family_asv[name]} ASVs, "
            f"{reads} reads "
            f"({100 * reads / total_reads:.2f}%)\n"
        )

    out.write("\n=== TOP 20 GENERA BY READS ===\n")

    for name, reads in genus_reads.most_common(20):
        out.write(
            f"{name}: {genus_asv[name]} ASVs, "
            f"{reads} reads "
            f"({100 * reads / total_reads:.2f}%)\n"
        )

    out.write("\n=== SPECIES ASSIGNMENTS ===\n")

    for name, reads in species_reads.most_common():
        if name != "Unassigned":
            out.write(
                f"{name}: {species_asv[name]} ASVs, "
                f"{reads} reads\n"
            )

    out.write("\n=== AMBIGUITY ===\n")
    out.write(
        f"Ambiguous genus: {len(ambiguous_genus_asvs)} ASVs, "
        f"{sum(int(x['Total_reads']) for x in ambiguous_genus_asvs)} reads\n"
    )

    out.write(
        f"Ambiguous species: {len(ambiguous_species_asvs)} ASVs, "
        f"{sum(int(x['Total_reads']) for x in ambiguous_species_asvs)} reads\n"
    )

print("Taxonomy + abundance QC completed.")
print("Joined table:", OUT_TABLE)
print("Summary:", OUT_SUMMARY)
