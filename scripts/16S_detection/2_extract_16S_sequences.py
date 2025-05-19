##!/usr/bin/env python

import os
import argparse
from Bio import SeqIO
from Bio.Seq import Seq

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-g', '--gff_root', required=True, help="Dossier GFF+FASTA généré par barrnap")
    parser.add_argument('-o', '--output_dir', required=True, help="Dossier de sortie pour séquences 16S extraites")
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    PADDING = 100  # contexte de 100pb autour

    for binner in os.listdir(args.gff_root):
        for sample in os.listdir(os.path.join(args.gff_root, binner)):
            sample_dir = os.path.join(args.gff_root, binner, sample)
            fasta_files = [f for f in os.listdir(sample_dir) if f.endswith('.fasta')]

            if not fasta_files:
                print(f"⚠️  FASTA  manquant pour {binner}/{sample}")
                continue

            output_fasta = os.path.join(args.output_dir, f"{binner}_{sample}_16S_extracted.fasta")
            with open(output_fasta, 'w') as out_f:
                for fasta_file in fasta_files:
                    fasta_path = os.path.join(sample_dir, fasta_file)
                    gff_path = fasta_path.replace('.fasta', '.gff')

                    if not os.path.exists(gff_path) or os.path.getsize(gff_path) == 0:
                        print(f"⚠️  GFF manquant pour {binner}/{sample}")
                        continue

                    records = {record.id: record for record in SeqIO.parse(fasta_path, 'fasta')}
                    with open(gff_path) as gf:
                        for line in gf:
                            if line.startswith('#') or line.strip() == '':
                                continue
                            fields = line.strip().split('\t')
                            if len(fields) < 9:
                                continue
                            if fields[2] == 'rRNA' and '16S' in fields[8]:
                                scaffold_id = fields[0]
                                start = int(fields[3])
                                end = int(fields[4])
                                strand = fields[6]

                                if scaffold_id not in records:
                                    continue

                                seq = records[scaffold_id].seq
                                extract_start = max(0, start - 1 - PADDING)
                                extract_end = min(len(seq), end + PADDING)
                                extracted_seq = seq[extract_start:extract_end]

                                if strand == '-':
                                    extracted_seq = extracted_seq.reverse_complement()

                                out_f.write(f">{scaffold_id}_{start}_{end}_{strand}\n{extracted_seq}\n")

    print("\n✅ Extraction 16S terminée.")

if __name__ == "__main__":
    main()
