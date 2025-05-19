#!/usr/bin/env python

# === RECALCUL SUMMARY_16S à partir des GFF existants ===

import os
import pandas as pd
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--gff_root', required=True, help="Dossier racine contenant {binner}/{sample}/*.gff")
parser.add_argument('--output_summary', required=True, help="Fichier TSV de sortie pour le nouveau summary")
args = parser.parse_args()

binner_samples = []
for binner in os.listdir(args.gff_root):
    binner_path = os.path.join(args.gff_root, binner)
    if not os.path.isdir(binner_path):
        continue

    for sample in os.listdir(binner_path):
        sample_path = os.path.join(binner_path, sample)
        if not os.path.isdir(sample_path):
            continue

        gff_files = [f for f in os.listdir(sample_path) if f.endswith('.gff')]

        n_total = len(gff_files)
        n_with_16s = 0

        for gff_file in gff_files:
            gff_path = os.path.join(sample_path, gff_file)
            with open(gff_path, 'r') as f:
                for line in f:
                    if not line.startswith('#') and '16S' in line:
                        n_with_16s += 1
                        break  # On compte un bin une seule fois

        percent = (n_with_16s / n_total * 100) if n_total else 0
        binner_samples.append({
            'binner': binner,
            'sample': sample,
            'n_bins_total': n_total,
            'n_bins_with_16s': n_with_16s,
            'percent_with_16s': round(percent, 2),
        })

df = pd.DataFrame(binner_samples)
df.to_csv(args.output_summary, sep='\t', index=False)

print(f"✅ Nouveau résumé sauvegardé dans {args.output_summary}")
