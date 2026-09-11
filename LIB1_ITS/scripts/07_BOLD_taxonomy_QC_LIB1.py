import csv
from collections import defaultdict, Counter

CLASS_FILE = "taxonomy/processed/LIB1_BOLD_classifications_all.csv"
HITS_FILE  = "taxonomy/processed/LIB1_BOLD_combined_hits_all.csv"

OUT_TABLE = "taxonomy/qc/LIB1_BOLD_taxonomy_QC.csv"
OUT_SUMMARY = "taxonomy/qc/LIB1_BOLD_taxonomy_QC_summary.txt"

CONF_THRESHOLD = 80.0

# ---------------------------------------------------------
# 1. Read combined hits
# ---------------------------------------------------------

hits_by_asv = defaultdict(list)

with open(HITS_FILE, newline="", encoding="utf-8-sig") as f:
    reader = csv.DictReader(f)

    for row in reader:
        asv = row["Query ID"]

        try:
            identity = float(row["ID%"])
        except:
            identity = None

        row["_identity"] = identity
        hits_by_asv[asv].append(row)

# ---------------------------------------------------------
# 2. Read classifications and build QC table
# ---------------------------------------------------------

results = []

with open(CLASS_FILE, newline="", encoding="utf-8-sig") as f:
    reader = csv.DictReader(f)

    for row in reader:

        asv = row["Query ID"]

        rank = row["Tax Rank"].strip()
        confidence_raw = row["Confidence"].strip()

        try:
            confidence = float(confidence_raw)
        except:
            confidence = None

        asv_hits = hits_by_asv.get(asv, [])

        valid_ids = [
            h["_identity"]
            for h in asv_hits
            if h["_identity"] is not None
        ]

        best_identity = max(valid_ids) if valid_ids else None

        # Hits tied at the maximum identity
        top_hits = [
            h for h in asv_hits
            if h["_identity"] == best_identity
        ] if best_identity is not None else []

        top_families = sorted(set(
            h["Family"] for h in top_hits if h["Family"].strip()
        ))

        top_genera = sorted(set(
            h["Genus"] for h in top_hits if h["Genus"].strip()
        ))

        top_species = sorted(set(
            h["Species"] for h in top_hits if h["Species"].strip()
        ))

        # Assignment status
        if not rank:
            assignment_status = "Unassigned"
        elif confidence is not None and confidence < CONF_THRESHOLD:
            assignment_status = "Low-confidence"
        else:
            assignment_status = "Assigned"

        # Ambiguity among equally best hits
        ambiguous_family = "YES" if len(top_families) > 1 else "NO"
        ambiguous_genus = "YES" if len(top_genera) > 1 else "NO"
        ambiguous_species = "YES" if len(top_species) > 1 else "NO"

        results.append({
            "ASV": asv,
            "Tax_Rank": rank if rank else "Unassigned",
            "Phylum": row["Phylum"],
            "Class": row["Class"],
            "Order": row["Order"],
            "Family": row["Family"],
            "Subfamily": row["Subfamily"],
            "Genus": row["Genus"],
            "Species": row["Species"],
            "Confidence": confidence_raw,
            "Supporting_Recs": row["Supporting Recs"],
            "Best_ID_percent": (
                f"{best_identity:.2f}"
                if best_identity is not None else ""
            ),
            "Number_of_BOLD_hits": len(asv_hits),
            "Top_hit_families": ";".join(top_families),
            "Top_hit_genera": ";".join(top_genera),
            "Top_hit_species": ";".join(top_species),
            "Ambiguous_family": ambiguous_family,
            "Ambiguous_genus": ambiguous_genus,
            "Ambiguous_species": ambiguous_species,
            "Assignment_status": assignment_status
        })

# ---------------------------------------------------------
# 3. Write QC table
# ---------------------------------------------------------

fieldnames = list(results[0].keys())

with open(OUT_TABLE, "w", newline="", encoding="utf-8") as f:
    writer = csv.DictWriter(f, fieldnames=fieldnames)
    writer.writeheader()
    writer.writerows(results)

# ---------------------------------------------------------
# 4. Generate summary
# ---------------------------------------------------------

rank_counts = Counter(r["Tax_Rank"] for r in results)
status_counts = Counter(r["Assignment_status"] for r in results)

family_counts = Counter(
    r["Family"] for r in results if r["Family"]
)

genus_counts = Counter(
    r["Genus"] for r in results if r["Genus"]
)

species_counts = Counter(
    r["Species"] for r in results if r["Species"]
)

amb_genus = sum(
    1 for r in results if r["Ambiguous_genus"] == "YES"
)

amb_species = sum(
    1 for r in results if r["Ambiguous_species"] == "YES"
)

with open(OUT_SUMMARY, "w", encoding="utf-8") as out:

    out.write("=== LIB1 BOLD TAXONOMY QC ===\n\n")

    out.write(f"Total ASVs: {len(results)}\n")
    out.write(f"Confidence threshold: {CONF_THRESHOLD}%\n\n")

    out.write("=== ASSIGNMENT STATUS ===\n")
    for k, v in status_counts.most_common():
        out.write(f"{k}: {v}\n")

    out.write("\n=== TAXONOMIC RANK ===\n")
    for k, v in rank_counts.most_common():
        out.write(f"{k}: {v}\n")

    out.write("\n=== AMBIGUITY AMONG BEST HITS ===\n")
    out.write(f"Ambiguous genus: {amb_genus}\n")
    out.write(f"Ambiguous species: {amb_species}\n")

    out.write("\n=== TOP 15 FAMILIES ===\n")
    for name, n in family_counts.most_common(15):
        out.write(f"{name}: {n}\n")

    out.write("\n=== TOP 15 GENERA ===\n")
    for name, n in genus_counts.most_common(15):
        out.write(f"{name}: {n}\n")

    out.write("\n=== SPECIES ASSIGNMENTS ===\n")
    for name, n in species_counts.most_common():
        out.write(f"{name}: {n}\n")

print("Taxonomy QC completed.")
print("QC table:", OUT_TABLE)
print("Summary:", OUT_SUMMARY)
