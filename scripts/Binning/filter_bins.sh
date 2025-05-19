#!/bin/bash
#SBATCH --job-name=filter_bin
#SBATCH --output=logs/filter_bin_%j.out
#SBATCH --error=logs/filter_bin_stats_%j.err
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G

echo "🔁 Activation de l’environnement Conda"
source activate mahabio_env

# Répertoires
SCRIPTS_DIR="/shared/home/asandri/MAHABIO_analysis/scripts/Binning"
STATS_DIR="/shared/home/asandri/MAHABIO_analysis/results/binning/stats"
MERGED_FILE="$STATS_DIR/all_binners_stats.tsv"
CHECKM_RESULTS="/shared/home/asandri/MAHABIO_analysis/results/binning/checkm"


python "$SCRIPTS_DIR/filter_bins.py" \
  -i "$STATS_DIR/all_binners_stats.tsv" \
  -b "$CHECKM_RESULTS" \
  -o "$STATS_DIR"

