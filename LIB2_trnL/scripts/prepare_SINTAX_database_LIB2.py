import csv

fasta_file = "taxonomy/SINTAX/reference/Taberlet_trnL/derep_and_clean_db/Taberlet_trnL_derep_and_clean.fasta"
taxonomy_file = "taxonomy/SINTAX/reference/Taberlet_trnL/derep_and_clean_db/Taberlet_trnL_derep_and_clean_taxonomy.txt"
output_file = "taxonomy/SINTAX/reference/Taberlet_trnL_SINTAX.fasta"

ranks = ["k", "p", "c", "o", "f", "g", "s"]

taxonomy = {}

with open(taxonomy_file, newline="") as f:
    reader = csv.reader(f, delimiter="\t", quotechar='"')
    for row in reader:
        if len(row) < 2:
            continue

        seq_id = row[0]
        taxa = row[1].split(";")

        formatted = []

        for rank, taxon in zip(ranks, taxa):
            taxon = taxon.strip()

            if taxon and taxon != "NA":
                # Avoid whitespace problems in SINTAX FASTA headers
                taxon = taxon.replace(" ", "_")
                formatted.append(f"{rank}:{taxon}")

        taxonomy[seq_id] = ",".join(formatted)

written = 0
missing = 0

with open(fasta_file) as inp, open(output_file, "w") as out:
    for line in inp:
        if line.startswith(">"):
            seq_id = line[1:].strip()

            if seq_id in taxonomy:
                out.write(f">{seq_id};tax={taxonomy[seq_id]};\n")
                written += 1
            else:
                out.write(line)
                missing += 1
        else:
            out.write(line)

print("Reference sequences written:", written)
print("Headers without taxonomy:", missing)
print("Output:", output_file)
