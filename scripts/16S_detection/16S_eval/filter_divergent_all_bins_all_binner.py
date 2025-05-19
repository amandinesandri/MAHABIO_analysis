#!/usr/bin/env python3

import os
import argparse
from Bio import SeqIO
from itertools import combinations
from difflib import SequenceMatcher

def sequence_identity(seq1, seq2):
    return SequenceMatcher(None, str(seq1), str(seq2)).ratio()

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--fasta_dir', required=True, help="Répertoire contenant les fichiers FASTA 16S par bin")
    parser.add_argument('--threshold', type=float, default=0.97, help="Seuil d'identité min (ex: 0.97)")
    args = parser.parse_args()

    flagged_bins = []

    for file in os.listdir(args.fasta_dir):
            if not file.endswith("_16S_extracted.fasta"):
            continue


        records = list(SeqIO.parse(os.path.join(args.fasta_dir, file), "fasta"))
        if len(records) <= 1:
            print(f"⚠️ record:  {records} ")
            continue  # rien à comparer

        for r1, r2 in combinations(records, 2):
            print(f"⚠️ r1 du record:  {r1} ")
            print(f"⚠️ r2 du record:  {r2} ")
            iden = sequence_identity(r1.seq, r2.seq)
            print(f"⚠️ calcul sequence_identity:  {iden} ")
            if iden < args.threshold:
                flagged_bins.append((file, r1.id, r2.id, iden))
                # Ne PAS faire break : on veut toutes les paires divergentes !

    if flagged_bins:
        print("\n⚠️ Bins avec séquences 16S divergentes :")
        for fb in flagged_bins:
            print(f"{fb[0]} → {fb[1]} vs {fb[2]} = {fb[3]*100:.2f}% identité")
    else:
        print("✅ Aucun bin avec divergence de 16S détectée.")

if __name__ == "__main__":
    main()