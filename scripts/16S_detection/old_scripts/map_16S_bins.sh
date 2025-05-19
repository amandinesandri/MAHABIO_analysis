#!/bin/bash
#SBATCH --job-name=map_16S_bins
#SBATCH --output=logs/map_16S_bins_%j.out
#SBATCH --error=logs/map_16S_bins_%j.err
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G

echo "🚀 Activation de l’environnement Conda"
source activate mahabio_env

# === Répertoires ===
SCRIPT="/shared/home/asandri/MAHABIO_analysis/scripts/Binning/map_16S_bins.py"
INPUT_DIR="/shared/home/asandri/MAHABIO_analysis/results/binning/stats/16S_detection/16S_extracted"
OUTPUT_FASTA="/shared/home/asandri/MAHABIO_analysis/results/binning/stats/16S_detection/bins_individuels_16S.fasta"
OUTPUT_MAPPING="/shared/home/asandri/MAHABIO_analysis/results/binning/stats/16S_detection/bins_16S_mapping.tsv"

mkdir -p logs

echo "🔎 Lancement du mapping bin -> 16S..."
python "$SCRIPT" \
  -i "$INPUT_DIR" \
  -o "$OUTPUT_FASTA" \
  -m "$OUTPUT_MAPPING"

echo "✅ Terminé : résultats dans $OUTPUT_FASTA et $OUTPUT_MAPPING"
