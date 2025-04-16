#!/bin/bash
#SBATCH --job-name=coverage
#SBATCH --output=logs/coverage_%A_%a.out
#SBATCH --error=logs/coverage_%A_%a.err
#SBATCH --time=24:00:00
#SBATCH --mem=4G
#SBATCH --array=0-2

module load samtools

# Définition des fichiers BAM
bam_files=(
  "/shared/home/asandri/MAHABIO_analysis/results/coverage/Sj_contigs_more_than_300bp_cleaned_sorted.bam"
  "/shared/home/asandri/MAHABIO_analysis/results/coverage/Fk_contigs_more_than_300bp_cleaned_sorted.bam"
  "/shared/home/asandri/MAHABIO_analysis/results/coverage/C_contigs_more_than_300bp_cleaned_sorted.bam"
)

# Sélection du fichier BAM en fonction de l'ID de la tâche SLURM
BAM="${bam_files[$SLURM_ARRAY_TASK_ID]}"
BASENAME=$(basename "$BAM" .bam)

echo "Traitement de $BAM..."

# Indexation si nécessaire
if [ ! -f "${BAM}.bai" ]; then
    samtools index "$BAM"
fi

# Calcul de la couverture
samtools depth "$BAM" | awk '{sum[$1]+=$3} END {for (scaffold in sum) print scaffold, sum[scaffold]}' > "${BASENAME}.depth.txt"
