#!/bin/bash
#SBATCH --job-name=extract_16S
#SBATCH --output=logs/extract_16S_%A_%a.out
#SBATCH --error=logs/extract_16S_%A_%a.err
#SBATCH --array=00  # sera remplacé dynamiquement par extract_all_16S.sh
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G

# ------------------ ⏱️ DÉBUT ------------------
START_TIME=$(date +%s)
echo "⏱️ Job lancé le : $(date '+%Y-%m-%d %H:%M:%S')"
# ------------------------------------------------

source activate mahabio_env

SAMPLES_LIST=$1
sample=$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" $SAMPLES_LIST)

binner=$(echo $sample | cut -d'/' -f1)
sample_name=$(echo $sample | cut -d'/' -f2)

echo "🔍 Extraction pour $binner/$sample_name"

python extract_16S_sequences.py \
    -g /shared/home/asandri/MAHABIO_analysis/results/binning/Affiliation/16S_detection \
    -o /shared/home/asandri/MAHABIO_analysis/results/binning/Affiliation/16S_extracted
# ------------------ ✅ FIN ------------------
END_TIME=$(date +%s)
DURATION=$((END_TIME - START_TIME))
echo "✅ Job terminé le : $(date '+%Y-%m-%d %H:%M:%S')"

hours=$((DURATION / 3600))
minutes=$(((DURATION % 3600) / 60))
seconds=$((DURATION % 60))
echo "⏳ Durée totale : ${hours}h ${minutes}min ${seconds}s"
# ------------------------------------------------