#!/bin/bash
#SBATCH --job-name=bin_stats_summary
#SBATCH --output=logs/bin_stats_summary_%j.out
#SBATCH --error=logs/bin_stats_summary_%j.err
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G

# === Activer Conda ===
echo "🔁 Activation de l’environnement conda"
source activate mahabio_env

# === Répertoires ===
SCRIPT="/shared/home/asandri/MAHABIO_analysis/scripts/Binning/bin_stats_summary.py"
RESULT_DIR="/shared/home/asandri/MAHABIO_analysis/results/binning/checkm"
OUT_DIR="/shared/home/asandri/MAHABIO_analysis/results/binning/stats"

mkdir -p "$OUT_DIR" logs

# === Binners à traiter ===
BINS=("semibin" "maxbin" "vamb")

# === Étape 1 : générer les fichiers de stats individuels ===
for binner in "${BINS[@]}"; do
    echo "📊 Traitement de $binner..."

    find "$RESULT_DIR/$binner" -mindepth 2 -type f -name "quality_summary.tsv" | while read file; do
        SAMPLE_DIR=$(dirname "$file")
        SAMPLE_NAME=$(basename "$SAMPLE_DIR")

        echo "Traitement de $binner - Sample détecté : $SAMPLE_NAME"

        python "$SCRIPT" \
            -i "$SAMPLE_DIR" \
            -o "$OUT_DIR/${binner}_${SAMPLE_NAME}_stats.tsv"
    done
done

