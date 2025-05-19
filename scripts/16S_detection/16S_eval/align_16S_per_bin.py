#!/usr/bin/env python3

import os
import subprocess
import argparse
from Bio import SeqIO

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--fasta_dir', required=True, help="Répertoire contenant les fichiers FASTA extraits")
    parser.add_argument('--output_dir', required=True, help="Répertoire de sortie des alignements MAFFT")
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    for root, _, files in os.walk(args.fasta_dir):
        for file in files:
            if not file.endswith("_16S_extracted.fasta"):
                continue
            fasta_path = os.path.join(root, file)

            # Vérifier nombre de séquences
            seqs = list(SeqIO.parse(fasta_path, "fasta"))
            if len(seqs) < 2:
                print(f"⚠️  Moins de 2 séquences dans {file}, alignement MAFFT ignoré.")
                continue

            bin_id = file.replace("_16S_extracted.fasta", "")
            output_path = os.path.join(args.output_dir, f"{bin_id}_aligned.fasta")

            print(f"🧬 Alignement MAFFT : {file}")
            cmd = f"mafft --auto {fasta_path} > {output_path}"
            subprocess.run(cmd, shell=True)
            print(f"✅ Fichier aligné sauvegardé dans : {output_path}")

if __name__ == "__main__":
    main()
