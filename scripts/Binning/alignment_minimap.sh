#!/bin/bash
#SBATCH --job-name=alignment_minimap
#SBATCH --output=logs/alignment_minimap_%A_%a.out
#SBATCH --error=logs/alignment_minimap_%A_%a.err
#SBATCH --time=24:00:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=8
#SBATCH --array=0-1

module load samtools
module load bamtools
module load bwa
module load minimap2

# Définir les tableaux de contigs et de reads
contigs=(
    "/shared/home/asandri/MAHABIO/data/Sj_contigs_more_than_300bp_cleaned.fasta"
    "/shared/home/asandri/MAHABIO/data/Fk_contigs_more_than_300bp_cleaned.fasta"
)
reads=(
    "/shared/home/asandri/MAHABIO/data/Sj.fastq.gz"
    "/shared/home/asandri/MAHABIO/data/Fk.fastq.gz"
)

# Sélectionner les fichiers en fonction de l'ID de tâche
CTG="${contigs[$SLURM_ARRAY_TASK_ID]}"
READS="${reads[$SLURM_ARRAY_TASK_ID]}"
BASENAME=$(basename "$CTG" .fasta)

# Créer et se déplacer dans le répertoire de travail
WORKDIR="/shared/home/asandri/MAHABIO_analysis/results/coverage"
mkdir -p "$WORKDIR"
cd "$WORKDIR" || exit 1

echo "### Alignement de $READS sur $CTG avec minimap2"
minimap2 -ax map-ont -t 8 "$CTG" "$READS" | samtools view -bS - > "${BASENAME}.bam"

echo "### Tri et indexation"
samtools sort "${BASENAME}.bam" -o "${BASENAME}_sorted.bam"
samtools index "${BASENAME}_sorted.bam"

# Nettoyage des fichiers intermédiaires
#rm "${BASENAME}.bam"
