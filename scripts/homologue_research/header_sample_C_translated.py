from Bio import SeqIO
import sys

input_faa = sys.argv[1]
output_faa = sys.argv[2]

unique_records = []
id_counts = {}

for record in SeqIO.parse(input_faa, "fasta"):
    base_id = record.id.split()[0]  # Garder uniquement C|ctg000000001
    id_counts[base_id] = id_counts.get(base_id, 0) + 1
    new_id = f"{base_id}_seq{id_counts[base_id]}"
    record.id = new_id
    record.description = ""  # Supprimer description qui peut gêner
    unique_records.append(record)

SeqIO.write(unique_records, output_faa, "fasta")
