#!/usr/bin/env python

# === SCRIPT : FUSION AVEC COLONNES BINNER + SAMPLE SEPAREES ===

import os
import glob
import pandas as pd
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input_dir', required=True, help='Répertoire contenant tous les fichiers *_stats.tsv')
parser.add_argument('-o', '--output_file', required=True, help='Fichier TSV de sortie fusionné')
args = parser.parse_args()

# === Chercher tous les fichiers *_stats.tsv ===
tsv_files = glob.glob(os.path.join(args.input_dir, '*_stats.tsv'))

# === Fusionner les fichiers ===
all_dfs = []
for file in tsv_files:
    df = pd.read_csv(file, sep='\t', engine='python')

    # Déduire binner et sample à partir du nom du fichier
    basename = os.path.basename(file).replace('_stats.tsv', '')
    parts = basename.split('_', 1)
    if len(parts) != 2:
        raise ValueError(f"Format de fichier inattendu : {basename}, impossible de déduire binner et sample.")

    binner, sample = parts

    # Ajouter deux colonnes séparées
    df.insert(0, 'sample', sample)
    df.insert(0, 'binner', binner)

    all_dfs.append(df)

# === Concaténer tous les DataFrames ===
merged_df = pd.concat(all_dfs, ignore_index=True)

# === Sauvegarde
merged_df.to_csv(args.output_file, sep='\t', index=False)

print(f"\U0001F389 Fichier fusionné exporté vers {args.output_file}")
