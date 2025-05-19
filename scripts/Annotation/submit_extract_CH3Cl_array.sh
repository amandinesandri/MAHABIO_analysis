#!/bin/bash
#SBATCH --job-name=extract_chloro_gene
#SBATCH --output=logs/extract_chloro_gene_%j.out
#SBATCH --error=logs/extract_chloro_gene_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=1:00:00

# Fichier TSV avec les bins sélectionnés
BIN_FILE="/shared/home/asandri/MAHABIO_analysis/results/MAG_mapping/Bins_filtered_over_500kb.tsv"

# Calcul automatique du nombre de lignes (en ignorant l'en-tête)
N_JOBS=$(($(wc -l < "$BIN_FILE") - 1))

# Soumission du job array
sbatch --array=0-$N_JOBS extract_CH3Cl_array_job.sh "$BIN_FILE"
