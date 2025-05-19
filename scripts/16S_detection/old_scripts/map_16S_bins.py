#!/usr/bin/env python

import os
import argparse
import hashlib
from Bio import SeqIO

def hash_seq(seq):
    return hashlib.md5(str(seq).encode()).hexdigest()

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input_dir', required=True, help='Dossier contenant les *.fasta extraits')
    parser.add_argument('-o', '--output_fasta', required=True, help='Fichier FASTA de sortie par bin')
    parser.add_argument('-m', '--output_mapping', required=True, help='Fichier TSV de mapping binner-sample-bin -> sequence')
    args = parser.parse_args()

    records = []
    mapping = []

    for fasta_file in os.listdir(args.input_dir):
        if not fasta_file.endswith('.fasta'):
            continue

        binner, sample, *_ = fasta_file.replace('_16S_extracted.fasta', '').split('_', 2)
        fasta_path = os.path.join(args.input_dir, fasta_file)

        for record in SeqIO.parse(fasta_path, "fasta"):
            header_parts = record.id.split('_')
            bin_id = header_parts[0]
            start = header_parts[1]
            end = header_parts[2]
            strand = header_parts[3]

            seq_hash = hash_seq(record.seq)

            mapping.append({
                "binner": binner,
                "sample": sample,
                "bin_id": bin_id,
                "start": start,
                "end": end,
                "strand": strand,
                "sequence_length": len(record.seq),
                "seq_hash": seq_hash
            })

            record.id = f"{binner}|{sample}|{bin_id}|{start}-{end}({strand})"
            record.description = ""
            records.append(record)

    # Sauvegarder FASTA
    SeqIO.write(records, args.output_fasta, "fasta")

    # Sauvegarder TSV
    import pandas as pd
    df = pd.DataFrame(mapping)
    df.to_csv(args.output_mapping, sep='\t', index=False)

    print(f"\n✅ Fichiers créés : {args.output_fasta} et {args.output_mapping}")

if __name__ == "__main__":
    main()
