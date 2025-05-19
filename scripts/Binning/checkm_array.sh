#!/bin/bash
#SBATCH --job-name=checkm_array
#SBATCH --output=logs/checkm_array_%A_%a.out
#SBATCH --error=logs/checkm_array_%A_%a.err
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --array=0-2

ENV_NAME="mahabio_env"
BASE_INPUT="/shared/home/asandri/MAHABIO_analysis/results/binning"
BASE_OUTPUT="/shared/home/asandri/MAHABIO_analysis/results/binning/checkm"
BINNERS=("maxbin" "semibin" "vamb")

# Activer l'environnement
source ~/miniforge3/etc/profile.d/conda.sh
conda activate "$ENV_NAME"

BIN_NAME=${BINNERS[$SLURM_ARRAY_TASK_ID]}
INPUT_DIR="${BASE_INPUT}/${BIN_NAME}"
OUTPUT_DIR="${BASE_OUTPUT}/${BIN_NAME}"

echo "➡️  Traitement de $BIN_NAME"

mkdir -p "$OUTPUT_DIR"
checkm lineage_wf -x fa --threads 8 "$INPUT_DIR" "$OUTPUT_DIR"
checkm qa --tab_table --out_format 2 "$OUTPUT_DIR/lineage.ms" "$OUTPUT_DIR" > "$OUTPUT_DIR/quality_summary.tsv"

echo "✅ CheckM terminé pour $BIN_NAME"
