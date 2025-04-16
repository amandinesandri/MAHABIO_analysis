#!/bin/bash
#SBATCH --job-name=alignment_bwa
#SBATCH --output=logs/alignment_bwa_%j.out
#SBATCH --error=logs/alignment_bwa_%j.err
#SBATCH --time=24:00:00
#SBATCH --mem=4G

module load samtools
module load bamtools
module load bwa
module load minimap2
# Liste des fichiers BAM

echo $PWD 
WORKDIR="/shared/home/asandri/MAHABIO_analysis/results/coverage"
mkdir -p "$WORKDIR"
cd "$WORKDIR"

# Liste des fichiers

# # Boucle sur les fichiers
CTG="/shared/home/asandri/MAHABIO/data/C_contigs_more_than_300bp_cleaned.fasta"
#CTG="/shared/home/asandri/MAHABIO/data/C_contigs_more_than_300bp.fasta"
R1="/shared/home/asandri/MAHABIO/data/C_R1.fastq.gz"
R2="/shared/home/asandri/MAHABIO/data/C_R2.fastq.gz"

BASENAME=$(basename "$CTG" .fasta)

echo "### Indexation de $CTG avec BWA"
bwa index "$CTG"

echo "### Alignement des reads de $BASENAME"
bwa mem -t 8 "$CTG" "$R1" "$R2" | samtools view -bS - > "${BASENAME}.bam"

echo "### Tri et indexation"
samtools sort "${BASENAME}.bam" -o "${BASENAME}_sorted.bam"
samtools index "${BASENAME}_sorted.bam"

# Nettoyage si besoin
#mv "${BASENAME}_sorted.bam" "${BASENAME}.bam"
#mv "${BASENAME}_sorted.bam.bai" "${BASENAME}.bam.bai"
#rm "${BASENAME}.bam"

