import csv
import re
from collections import Counter

inp = "taxonomy/BOLD_SINTAX/results/LIB3_SINTAX_all.tsv"
out = "taxonomy/BOLD_SINTAX/processed/LIB3_SINTAX_taxonomy.csv"
summary_out = "taxonomy/BOLD_SINTAX/qc/LIB3_SINTAX_taxonomy_summary.txt"

rank_map = {
    "k": "Kingdom",
    "p": "Phylum",
    "c": "Class",
    "o": "Order",
    "f": "Family",
    "g": "Genus",
    "s": "Species"
}

rank_order = ["k", "p", "c", "o", "f", "g", "s"]

fields = ["ASV_ID"]

for r in rank_order:
    fields.append(rank_map[r])
    fields.append(rank_map[r] + "_confidence")

fields += [
    "Strand",
    "Deepest_rank",
    "Deepest_taxon",
    "Accepted_taxonomy"
]

summary = Counter()
written = 0
invalid = 0

with open(inp, encoding="utf-8") as fin, \
     open(out, "w", newline="", encoding="utf-8") as fout:

    writer = csv.DictWriter(fout, fieldnames=fields)
    writer.writeheader()

    for line in fin:
        line = line.rstrip("\n")

        if not line:
            continue

        parts = line.split("\t")

        if len(parts) < 4:
            invalid += 1
            continue

        asv_id = parts[0]
        raw_assignment = parts[1]
        strand = parts[2]
        accepted = parts[3]

        row = {field: "" for field in fields}
        row["ASV_ID"] = asv_id
        row["Strand"] = strand
        row["Accepted_taxonomy"] = accepted

        # Parse taxonomy + confidence from SINTAX raw assignment
        matches = re.findall(
            r'([kpcofgs]):([^,()]+)\(([0-9.]+)\)',
            raw_assignment
        )

        confidence_lookup = {}

        for prefix, taxon, confidence in matches:
            confidence_lookup[prefix] = (taxon, confidence)

        # Parse taxonomy that passed the 0.80 cutoff
        accepted_lookup = {}

        if accepted:
            for item in accepted.split(","):
                if ":" not in item:
                    continue

                prefix, taxon = item.split(":", 1)

                if prefix in rank_map:
                    accepted_lookup[prefix] = taxon

        deepest_prefix = None

        for prefix in rank_order:
            rank_name = rank_map[prefix]

            if prefix in accepted_lookup:
                taxon = accepted_lookup[prefix]

                row[rank_name] = taxon

                if prefix in confidence_lookup:
                    row[rank_name + "_confidence"] = confidence_lookup[prefix][1]

                deepest_prefix = prefix

        if deepest_prefix:
            row["Deepest_rank"] = rank_map[deepest_prefix]
            row["Deepest_taxon"] = accepted_lookup[deepest_prefix]

            summary[rank_map[deepest_prefix]] += 1
        else:
            row["Deepest_rank"] = "Unassigned"
            summary["Unassigned"] += 1

        writer.writerow(row)
        written += 1

with open(summary_out, "w", encoding="utf-8") as f:
    f.write("=== LIB3 SINTAX TAXONOMY SUMMARY ===\n")
    f.write("Confidence cutoff: 0.80\n")
    f.write(f"ASVs processed: {written}\n")
    f.write(f"Invalid input rows: {invalid}\n\n")

    f.write("Deepest accepted taxonomic rank:\n")

    for rank in [
        "Species",
        "Genus",
        "Family",
        "Order",
        "Class",
        "Phylum",
        "Kingdom",
        "Unassigned"
    ]:
        n = summary[rank]
        pct = 100 * n / written if written else 0
        f.write(f"{rank}: {n} ({pct:.2f}%)\n")

print("=== LIB3 SINTAX TAXONOMY TABLE COMPLETED ===")
print("ASVs processed:", written)
print("Invalid rows:", invalid)
print("Output:", out)
print("Summary:", summary_out)
