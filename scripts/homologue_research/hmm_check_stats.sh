#!/bin/bash
#SBATCH --job-name=hmm_stats
#SBATCH --output=logs/hmm_stats_%j.out
#SBATCH --error=logs/hmm_stats_%j.err
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G

# === 1. Activation de l'environnement Conda ===
source activate mahabio_env

# === 2. Définir le chemin du fichier combiné ===
COMBINED_TSV="/shared/home/asandri/MAHABIO_analysis/results/hmm/summary/combined_hits_summary.tsv"

# === 3. Lancer l'analyse avec un script Python inline ===
python - << EOF
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Charger les données
df = pd.read_csv("${COMBINED_TSV}", sep="\t", on_bad_lines='skip')
# AFFICHER LES NOMS DE COLONNES VRAIS
print("\n==== Colonnes détectées dans le fichier ====\n")
print(df.columns)
# Résumé statistique
print("\\n===== Résumé statistique des scores, evalues et coverage =====\\n")
print(df[['evalue', 'score', 'coverage']].describe())

# Histogramme des scores
plt.figure(figsize=(6,4))
plt.hist(df['score'], bins=50, color='skyblue')
plt.title('Distribution des scores HMM')
plt.xlabel('Score')
plt.ylabel('Nombre de hits')
plt.grid(True)
plt.tight_layout()
plt.savefig("/shared/home/asandri/MAHABIO_analysis/results/hmm/plots/score_distribution.png")

# Histogramme des e-values (log10 transformées)
plt.figure(figsize=(6,4))
plt.hist(np.log10(df['evalue'].replace(0, 1e-300)), bins=50, color='salmon')  # éviter log(0)
plt.title('Distribution des e-values (log10)')
plt.xlabel('log10(e-value)')
plt.ylabel('Nombre de hits')
plt.grid(True)
plt.tight_layout()
plt.savefig("/shared/home/asandri/MAHABIO_analysis/results/hmm/plots/evalue_log10_distribution.png")

# Histogramme de la couverture
plt.figure(figsize=(6,4))
plt.hist(df['coverage'], bins=50, color='lightgreen')
plt.title('Distribution de la couverture')
plt.xlabel('Coverage')
plt.ylabel('Nombre de hits')
plt.grid(True)
plt.tight_layout()
plt.savefig("/shared/home/asandri/MAHABIO_analysis/results/hmm/plots/coverage_distribution.png")

print("\\n📊 Plots sauvegardés dans /results/hmm/plots/")
EOF

echo "✅ Analyse terminée. Vérifie le fichier .out pour le résumé et regarde les plots générés."
