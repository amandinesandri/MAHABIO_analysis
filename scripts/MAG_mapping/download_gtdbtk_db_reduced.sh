#!/bin/bash
#SBATCH --job-name=download_gtdbtk_reduced
#SBATCH --output=logs/download_gtdbtk_reduced.log
#SBATCH --error=logs/download_gtdbtk_reduce.err
#SBATCH --time=3-00:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=4G

# Active Conda et l'environnement contenant aria2
source /shared/home/asandri/miniforge3/etc/profile.d/conda.sh
source activate aria_env

# Dossier de destination
DEST_DIR="/shared/home/asandri/MAHABIO_analysis/ressources/gtdbtk_data/gtdbtk_msa_reduced"
mkdir -p "$DEST_DIR"
cd "$DEST_DIR" || exit 1

# URL de la base GTDB-Tk
URL="https://data.gtdb.ecogenomic.org/releases/release226/226.0/auxillary_files/gtdbtk_r226_msa_reduced.tar.gz"

echo "🚀 Démarrage du téléchargement GTDB-Tk reduced avec aria2..."
aria2c -x 8 -s 8 -c "$URL"

echo "📦 Extraction de l’archive..."
tar -xvzf gtdbtk_r226_msa_reduced.tar.gz -C gtdbtk_msa_reduced

echo "✅ Téléchargement terminé. Tu peux exporter GTDBTK_DATA_PATH :"
echo "export GTDBTK_DATA_reduced_PATH=$DEST_DIR/gtdbtk_msa_reduced"
