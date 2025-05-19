#!/usr/bin/env python

import os
import argparse
import matplotlib.pyplot as plt
from Bio import AlignIO

def plot_alignment(input_dir, output_dir):
    os.makedirs(output_dir, exist_ok=True)
    for aligned_file in os.listdir(input_dir):
        if not aligned_file.endswith("_aligned.fasta"):
            continue
        input_path = os.path.join(input_dir, aligned_file)
        try:
            alignment = AlignIO.read(input_path, "fasta")
        except Exception as e:
            print(f"⚠️ Erreur lecture {aligned_file} : {e}")
            continue

        plt.figure(figsize=(14, 0.4 * len(alignment)))
        y_labels = []
        for idx, record in enumerate(alignment):
            sequence = str(record.seq)
            colors = ["black" if base != "-" else "white" for base in sequence]
            plt.scatter(range(len(sequence)), [idx]*len(sequence), c=colors, marker='s', s=10)
            y_labels.append(record.id)

        plt.yticks(range(len(alignment)), y_labels, fontsize=8)
        plt.title(aligned_file.replace("_aligned.fasta", ""))
        plt.xlabel("Position dans l'alignement")
        plt.ylabel("Séquences (bin)")
        plt.tight_layout()
        output_plot = os.path.join(output_dir, aligned_file.replace("_aligned.fasta", "_alignment_tagged.png"))
        plt.savefig(output_plot, dpi=300)
        plt.close()
        print(f"✅ Visualisation sauvegardée : {output_plot}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input_dir', required=True)
    parser.add_argument('-o', '--output_dir', required=True)
    args = parser.parse_args()
    plot_alignment(args.input_dir, args.output_dir)
