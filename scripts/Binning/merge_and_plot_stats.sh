#!/bin/bash
#SBATCH --job-name=merge_plot_stats
#SBATCH --output=logs/merge_plot_stats_%j.out
#SBATCH --error=logs/merge_plot_stats_%j.err
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G

echo "🔁 Activation de l’environnement Conda"
source activate mahabio_env

# Répertoires
SCRIPTS_DIR="/shared/home/asandri/MAHABIO_analysis/scripts/Binning"
STATS_DIR="/shared/home/asandri/MAHABIO_analysis/results/binning/stats"
MERGED_FILE="$STATS_DIR/all_binners_stats.tsv"
PLOT_FILE="$STATS_DIR/checkm_quality_plot.png"

mkdir -p logs

# echo "📎 Fusion des fichiers *_stats.tsv..."
# python "$SCRIPTS_DIR/merge_bin_stats.py" \
#   -i "$STATS_DIR" \
#   -o "$MERGED_FILE"


echo "📊 Génération du graphique qualité CheckM..."
python "$SCRIPTS_DIR/plot_bin_stats.py"
# python "$SCRIPTS_DIR/plot_bin_stats.py" \
#   -i "$MERGED_FILE" \
#   -o "$PLOT_FILE"

echo "✅ Terminé : graphique dans $PLOT_FILE"
