#!/bin/bash
#SBATCH --job-name=16S_hits_summary
#SBATCH --output=logs/16S_hits_summary_%x_%j.out
#SBATCH --error=logs/16S_hits_summary_%x_%j.err
#SBATCH --time=3-00:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=4G

# ------------------ ⏱️ DÉBUT ------------------
START_TIME=$(date +%s)
echo "⏱️ Job lancé le : $(date '+%Y-%m-%d %H:%M:%S')"
# ------------------------------------------------

# (Optionnel) Activation de ton environnement Conda
source /shared/home/asandri//miniforge3/etc/profile.d/conda.sh
source activate mahabio_env

#python generate_16S_summary.py
python detect_shared_contigs_16S_hits.py

# ------------------ ✅ FIN ------------------
END_TIME=$(date +%s)
DURATION=$((END_TIME - START_TIME))
echo "✅ Job terminé le : $(date '+%Y-%m-%d %H:%M:%S')"

hours=$((DURATION / 3600))
minutes=$(((DURATION % 3600) / 60))
seconds=$((DURATION % 60))
echo "⏳ Durée totale : ${hours}h ${minutes}min ${seconds}s"
# ------------------------------------------------
