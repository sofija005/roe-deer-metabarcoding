import csv
import re
from collections import Counter

inp = "taxonomy/SINTAX/raw_results/LIB2_SINTAX_raw.tsv"
out = "taxonomy/SINTAX/processed/LIB2_SINTAX_taxonomy_QC.csv"
summary = "taxonomy/SINTAX/qc/LIB2_SINTAX_taxonomy_QC_summary.txt"

ranks = [
    ("k", "Kingdom"),
    ("p", "Phylum"),
    ("c", "Class"),
    ("o", "Order"),
    ("f", "Family"),
    ("g", "Genus"),
    ("s", "Species"),
]

def parse_raw_taxonomy(s):
    result = {}
    for part in s.split(","):
        m = re.match(r"([kpcofgs]):(.+)\(([\d.]+)\)$", part)
        if m:
            code, taxon, conf = m.groups()
            result[code] = (taxon.replace("_", " "), float(conf))
    return result

rows = []

with open(inp) as f:
    reader = csv.reader(f, delimiter="\t")

    for row in reader:
        asv = row[0]
        raw_tax = row[1] if len(row) > 1 else ""
        strand = row[2] if len(row) > 2 else ""
        accepted = row[3] if len(row) > 3 else ""

        parsed = parse_raw_taxonomy(raw_tax)

        accepted_codes = set()
        if accepted:
            for x in accepted.split(","):
                if ":" in x:
                    accepted_codes.add(x.split(":", 1)[0])

        rec = {
            "ASV": asv,
            "Strand": strand,
            "Raw_SINTAX": raw_tax,
            "Accepted_SINTAX": accepted,
        }

        accepted_rank = "Unassigned"
        accepted_conf = ""

        for code, name in ranks:
            taxon = ""
            conf = ""

            if code in parsed:
                raw_taxon, raw_conf = parsed[code]
                conf = raw_conf

                if code in accepted_codes:
                    taxon = raw_taxon
                    accepted_rank = name
                    accepted_conf = raw_conf

            rec[name] = taxon
            rec[name + "_confidence"] = conf

        rec["Accepted_rank"] = accepted_rank
        rec["Accepted_confidence"] = accepted_conf

        rows.append(rec)

fields = [
    "ASV", "Strand",
    "Kingdom", "Kingdom_confidence",
    "Phylum", "Phylum_confidence",
    "Class", "Class_confidence",
    "Order", "Order_confidence",
    "Family", "Family_confidence",
    "Genus", "Genus_confidence",
    "Species", "Species_confidence",
    "Accepted_rank", "Accepted_confidence",
    "Raw_SINTAX", "Accepted_SINTAX"
]

with open(out, "w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=fields)
    writer.writeheader()
    writer.writerows(rows)

rank_counts = Counter(r["Accepted_rank"] for r in rows)

with open(summary, "w") as f:
    f.write("LIB2 SINTAX TAXONOMY QC SUMMARY\n")
    f.write("===============================\n")
    f.write(f"Total ASVs: {len(rows)}\n")
    f.write("SINTAX confidence threshold: 0.80\n\n")

    f.write("Lowest accepted taxonomic rank:\n")
    for rank, n in rank_counts.most_common():
        f.write(f"{rank}: {n}\n")

print("Taxonomy table:", out)
print("Summary:", summary)
print("Total ASVs:", len(rows))
print("\nAccepted ranks:")
for rank, n in rank_counts.most_common():
    print(rank, n)
