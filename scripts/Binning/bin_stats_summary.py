#!/usr/bin/env python

import os
import pandas as pd
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input', required=True, help='Répertoire contenant quality_summary.tsv')
parser.add_argument('-o', '--output', required=True, help='Fichier de sortie TSV combiné')
args = parser.parse_args()

# === Définir le chemin
quality_path = os.path.join(args.input, 'quality_summary.tsv')

# === Vérification existence
if not os.path.exists(quality_path):
    raise FileNotFoundError(f"Fichier introuvable : {quality_path}")

# === Préparer à lire : ignorer les lignes INFO
valid_lines = []
with open(quality_path, 'r') as f:
    for line in f:
        if not line.startswith('['):  # ignorer les lignes de log [INFO]
            valid_lines.append(line)

# === Sauvegarder dans un fichier temporaire nettoyé
cleaned_path = os.path.join(args.input, 'quality_summary_cleaned.tsv')
with open(cleaned_path, 'w') as f:
    f.writelines(valid_lines)

# === Lire avec pandas en séparateurs flexibles (espaces multiples)
df = pd.read_csv(cleaned_path, sep='\t')
# === Vérifications colonnes
expected_cols = ['Bin Id', 'Completeness', 'Contamination']
for col in expected_cols:
    if col not in df.columns:
        raise ValueError(f"Colonne {col} manquante dans quality_summary.tsv")

# === Tri par Complétude décroissante
df = df.sort_values(by='Completeness', ascending=False)

# === Sauvegarde
df.to_csv(args.output, sep='\t', index=False)
print(f"✅ Résumé CheckM exporté vers {args.output}")
