#!/usr/bin/env python3

import os
import argparse
from Bio import SeqIO
from Bio import pairwise2
import pandas as pd

def compute_identity(seq1, seq2):
    alignments = pairwise2.align.globalxx(seq1, seq2, one_alignment_only=True, score_only=True)
    print(f"⚠️ Calcul alignments :   {alignments} ")
    max_len = max(len(seq1), len(seq2))
    print(f"⚠️ longueur max:  {max} ")
    return (alignments / max_len) * 100 if max_len > 0 else 0

def main(input_dir, output_dir):
    os.makedirs(output_dir, exist_ok=True)
    for aligned_file in os.listdir(input_dir):
        print(f"⚠️ aligned file utilisé:  {aligned_file} ")
        if not aligned_file.endswith("_aligned.fasta"):
            continue
        input_path = os.path.join(input_dir, aligned_file)
        print(f"⚠️ input path utilisé:  {input_path} ")
        records = list(SeqIO.parse(input_path, "fasta"))
        print(f"⚠️ le record ça donne quoi : :  {aligned_file} ")
        if len(records) < 2:
            print(f"⚠️ Pas assez de séquences dans {aligned_file} pour calculer une matrice.")
            continue

        matrix = []
        labels = [r.id for r in records]
        for r1 in records:
            print(f"⚠️ r1 du record:  {r1} ")
            row = []
            for r2 in records:
                print(f"⚠️ r2 du record:  {r2} ")
                identity = compute_identity(str(r1.seq), str(r2.seq))
                print(f"⚠️ Identity calcul :  {identity} ")
                row.append(round(identity, 2))
            matrix.append(row)

        df = pd.DataFrame(matrix, columns=labels, index=labels)
        output_matrix = os.path.join(output_dir, aligned_file.replace("_aligned.fasta", "_identity_matrix.tsv"))
        print(f"⚠️ output matrix pourquoi j' ai qu'un file ?:  {output_matrix} ")
        df.to_csv(output_matrix, sep="\t")
        print(f"✅ Matrice d'identité sauvegardée : {output_matrix}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input_dir', required=True)
    parser.add_argument('-o', '--output_dir', required=True)
    args = parser.parse_args()
    main(args.input_dir, args.output_dir)


