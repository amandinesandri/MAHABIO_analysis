#!/bin/bash
#SBATCH --job-name=align_16S
#SBATCH --output=logs/align_16S_%j.out
#SBATCH --error=logs/align_16S_%j.err
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G

echo "🚀 Activation de l’environnement Conda"
source activate mahabio_env

# Répertoires
INPUT_DIR="/shared/home/asandri/MAHABIO_analysis/results/binning/stats/16S_detection/16S_extracted"
OUTPUT_DIR="/shared/home/asandri/MAHABIO_analysis/results/binning/stats/16S_detection/16S_aligned"

mkdir -p "$OUTPUT_DIR"
mkdir -p logs

echo "📂 Lancement des alignements MAFFT..."
for fasta in "$INPUT_DIR"/*.fasta; do
    if [ -s "$fasta" ]; then
        base=$(basename "$fasta" .fasta)
        mafft --auto --thread 4 "$fasta" > "$OUTPUT_DIR/${base}_aligned.fasta"
        echo "✅ Alignement terminé pour $base"
    else
        echo "⚠️ Fichier vide ignoré : $fasta"
    fi
done

echo "🎯 Tous les alignements sont terminés."
