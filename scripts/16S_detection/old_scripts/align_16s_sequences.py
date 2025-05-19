#!/usr/bin/env python

# === SCRIPT 2 : ALIGNEMENT DES SEQUENCES 16S ===

import os
import argparse
import subprocess

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input_dir', required=True, help='Dossier contenant les fichiers FASTA extraits')
    parser.add_argument('-o', '--output_dir', required=True, help='Dossier de sortie pour les alignements MAFFT')
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    for file in os.listdir(args.input_dir):
        if not file.endswith("_16S_extracted.fasta"):
            continue

        input_fasta = os.path.join(args.input_dir, file)
        output_aln = os.path.join(args.output_dir, file.replace("_16S_extracted.fasta", "_aligned.fasta"))

        print(f"🔵 Alignement de {file}...")

        cmd = f"mafft --auto {input_fasta} > {output_aln}"
        result = subprocess.run(cmd, shell=True)

        if result.returncode != 0:
            print(f"⚠️  Erreur MAFFT sur {file}")

    print("\n✅ Alignements terminés pour toutes les séquences 16S !")

if __name__ == "__main__":
    main()
