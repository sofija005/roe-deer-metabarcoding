#!/usr/bin/env python3

import csv
import re
import argparse


def clean_taxon(value):
    """
    Clean taxonomic names so they can safely be used
    inside VSEARCH SINTAX FASTA headers.
    """
    value = value.strip()
    value = re.sub(r"[;,]", "_", value)
    value = re.sub(r"\s+", "_", value)
    return value


def main():
    parser = argparse.ArgumentParser(
        description="Prepare BOLD fungal ITS reference sequences for VSEARCH SINTAX."
    )

    parser.add_argument(
        "--input",
        required=True,
        help="BOLD full TSV export"
    )

    parser.add_argument(
        "--output",
        required=True,
        help="Output FASTA formatted for VSEARCH SINTAX"
    )

    args = parser.parse_args()

    ranks = [
        ("kingdom", "k"),
        ("phylum", "p"),
        ("class", "c"),
        ("order", "o"),
        ("family", "f"),
        ("genus", "g"),
        ("species", "s")
    ]

    accepted_markers = {"ITS", "ITS1", "ITS2"}

    total_records = 0
    marker_records = 0
    records_with_sequence = 0
    written = 0
    skipped_no_taxonomy = 0

    with open(
        args.input,
        encoding="utf-8",
        errors="replace"
    ) as fin, open(
        args.output,
        "w",
        encoding="utf-8"
    ) as fout:

        reader = csv.DictReader(fin, delimiter="\t")

        for row in reader:
            total_records += 1

            marker = row["marker_code"].strip()

            if marker not in accepted_markers:
                continue

            marker_records += 1

            seq = (
                row["nuc"]
                .strip()
                .upper()
                .replace(" ", "")
                .replace("-", "")
            )

            if not seq:
                continue

            records_with_sequence += 1

            taxonomy = []

            for column, prefix in ranks:
                value = clean_taxon(row[column])

                if value:
                    taxonomy.append(
                        f"{prefix}:{value}"
                    )

            if not taxonomy:
                skipped_no_taxonomy += 1
                continue

            process_id = row["processid"].strip()

            if not process_id:
                continue

            header = (
                f">{process_id};"
                f"tax={','.join(taxonomy)};"
            )

            fout.write(header + "\n")
            fout.write(seq + "\n")

            written += 1

    print("=== BOLD FUNGI ITS SINTAX DATABASE PREPARATION ===")
    print("Input records:", total_records)
    print("ITS/ITS1/ITS2 records:", marker_records)
    print("ITS records with sequence:", records_with_sequence)
    print("Skipped - no taxonomy:", skipped_no_taxonomy)
    print("SINTAX reference sequences written:", written)
    print("Output:", args.output)


if __name__ == "__main__":
    main()
