#!/bin/bash
#SBATCH --job-name=parse_cmuA_hits_param
#SBATCH --output=logs/parse_cmuA_hits_param_%j.out
#SBATCH --error=logs/parse_cmuA_hits_param_%j.err
#SBATCH --time=00:10:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1

# Activer l'environnement Conda
source activate mahabio_env

# Variables
INPUT_DIR="/shared/home/asandri/MAHABIO_analysis/results/binning/cmuA_hmmsearch"

# Paramètres de filtrage
EVAL_CUTOFF="1e-5"
SCORE_CUTOFF="50"

# Construire dynamiquement le nom du fichier de sortie
OUTPUT_FILE="/shared/home/asandri/MAHABIO_analysis/results/binning/cmuA_hmmsearch/cmuA_hits_summary_${EVAL_CUTOFF}_${SCORE_CUTOFF}.tsv"

# Lancer le parsing avec filtrage
python /shared/home/asandri/MAHABIO_analysis/scripts/homologue_research/binning/parse_cmuA_hits.py "$INPUT_DIR" "$OUTPUT_FILE" "$EVAL_CUTOFF" "$SCORE_CUTOFF"

echo "✅ Parsing terminé : résultats sauvegardés dans $OUTPUT_FILE"
