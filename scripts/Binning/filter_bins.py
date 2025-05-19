#!/usr/bin/env python

# === SCRIPT : DETECTION 16S PAR BINNER-SAMPLE AVEC AFFICHAGE DES CHEMINS ===

import os
import pandas as pd
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input', required=True, help='Fichier all_binners_stats.tsv')
parser.add_argument('-b', '--binning_root', required=True, help='Racine du dossier contenant les bins (ex: /shared/home/asandri/MAHABIO_analysis/results/binning/checkm)')
parser.add_argument('-o', '--output_dir', required=True, help='Dossier de sortie pour les résultats')
args = parser.parse_args()

os.makedirs(args.output_dir, exist_ok=True)

df = pd.read_csv(args.input, sep='\t')

if 'binner' not in df.columns or 'sample' not in df.columns:
    raise ValueError("Les colonnes 'binner' et 'sample' doivent être présentes dans le fichier d'entrée.")

binner_sample_list = df[['binner', 'sample']].drop_duplicates()

summary = []

def detect_16s_in_bin(bin_dir):
    has_16s = False
    gff_path = os.path.join(bin_dir, 'genes.gff')
    print(f"🔎 Chemin analysé : {gff_path}")
    if os.path.exists(gff_path):
        with open(gff_path) as gff_file:
            for line in gff_file:
                if '16S' in line:
                    has_16s = True
                    break
    else:
        print(f"⚠️ Fichier absent : {gff_path}")
    return has_16s

total_bins = 0
bins_with_16s = 0

for _, row in binner_sample_list.iterrows():
    binner = row['binner']
    sample = row['sample']
    binner_sample = f"{binner}_{sample}"
    bins_root = os.path.join(args.binning_root, binner, sample, 'bins')

    if not os.path.exists(bins_root):
        print(f"⚠️ Dossier bins introuvable pour {binner_sample} : {bins_root}")
        continue

    n_bins_total = 0
    n_bins_with_16s = 0
    detailed = []

    for bin_name in os.listdir(bins_root):
        bin_dir = os.path.join(bins_root, bin_name)
        if os.path.isdir(bin_dir):
            n_bins_total += 1
            total_bins += 1
            has_16s = detect_16s_in_bin(bin_dir)
            if has_16s:
                n_bins_with_16s += 1
                bins_with_16s += 1
            detailed.append({'bin': bin_name, 'has_16s': has_16s})

    percent_16s = (n_bins_with_16s / n_bins_total * 100) if n_bins_total > 0 else 0

    summary.append({
        'binner_sample': binner_sample,
        'n_bins_total': n_bins_total,
        'n_bins_with_16s': n_bins_with_16s,
        'percent_with_16s': round(percent_16s, 2)
    })

    pd.DataFrame(detailed).to_csv(
        os.path.join(args.output_dir, f'{binner_sample}_bins_16s_detail.tsv'),
        sep='\t', index=False
    )

pd.DataFrame(summary).to_csv(
    os.path.join(args.output_dir, 'bins_16s_summary_per_sample.tsv'),
    sep='\t', index=False
)

print("\n🎯 Analyse 16S terminée pour tous les binner-sample.")
print(f"\u2705 Bins analysés : {total_bins}")
print(f"\u2705 Bins avec 16S trouvés : {bins_with_16s}")
print(f"\u2705 Fichier résumé dans {os.path.join(args.output_dir, 'bins_16s_summary_per_sample.tsv')}")
