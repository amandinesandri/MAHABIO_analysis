#!/bin/bash
#SBATCH --job-name=download_gtdbtk
#SBATCH --output=logs/download_gtdbtk.log
#SBATCH --error=logs/download_gtdbtk.err
#SBATCH --time=3-00:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=4G

# Active Conda et l'environnement contenant aria2
source /shared/home/asandri/miniforge3/etc/profile.d/conda.sh
source activate aria_env

# Dossier de destination
DEST_DIR="/shared/home/asandri/MAHABIO_analysis/ressources/gtdbtk_data/release226"
mkdir -p "$DEST_DIR"
cd "$DEST_DIR" || exit 1

# URL de la base GTDB-Tk
URL="https://data.gtdb.ecogenomic.org/releases/release226/226.0/auxillary_files/gtdbtk_package/full_package/gtdbtk_r226_data.tar.gz"

echo "🚀 Démarrage du téléchargement GTDB-Tk R226 avec aria2..."
aria2c -x 8 -s 8 -c "$URL"

echo "📦 Extraction de l’archive..."
tar -xvzf gtdbtk_r226_data.tar.gz

echo "✅ Téléchargement terminé. Tu peux exporter GTDBTK_DATA_PATH :"
echo "export GTDBTK_DATA_PATH=$DEST_DIR/gtdbtk_r226_data"
