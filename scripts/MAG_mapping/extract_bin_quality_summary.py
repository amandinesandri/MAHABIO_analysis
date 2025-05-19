#!/usr/bin/env python3

import pandas as pd
import os

# === Fichiers d'entrée ===
stats_file = "/shared/home/asandri/MAHABIO_analysis/results/binning/stats/all_binners_stats.tsv"
filtered_bins_file = "/shared/home/asandri/MAHABIO_analysis/results/MAG_mapping/Bins_filtered_over_500kb.tsv"

# === Chargement des données ===
stats_df = pd.read_csv(stats_file, sep="\t")
filtered_df = pd.read_csv(filtered_bins_file, sep="\t")

# === Fusion des fichiers
merged = pd.merge(
    filtered_df,
    stats_df,
    how="left",
    left_on=["binner", "Sample", "filename"],
    right_on=["binner", "sample", "Bin Id"]
)

# === Ajout d’un label de qualité
def assess_quality(row):
    if row["Completeness"] >= 90 and row["Contamination"] < 5 and row["GC std (scaffolds > 1kbp)"] < 0.02:
        return "✅ Haute qualité"
    elif row["Completeness"] >= 70 and row["Contamination"] < 10 and row["GC std (scaffolds > 1kbp)"] < 0.03:
        return "⚠️ Moyenne qualité"
    else:
        return "❌ Faible qualité"

merged["Qualité"] = merged.apply(assess_quality, axis=1)

# === Colonnes finales à conserver
cols = [
    "binner", "Sample", "filename",
    "# genomes", "Completeness", "Contamination", "Strain heterogeneity",
    "Genome size (bp)", "# contigs", "N50 (scaffolds)", "N50 (contigs)",
    "Mean contig length (bp)", "Longest contig (bp)",
    "GC", "GC std (scaffolds > 1kbp)", "Coding density", "# predicted genes",
    "Qualité"
]

summary_df = merged[cols]

# === Enregistrement
outdir = "/shared/home/asandri/MAHABIO_analysis/results/MAG_mapping"
os.makedirs(outdir, exist_ok=True)
outfile = os.path.join(outdir, "Bins_filtered_quality_summary.tsv")
summary_df.to_csv(outfile, sep="\t", index=False)
print(f"✔️ Fichier généré : {outfile}")
