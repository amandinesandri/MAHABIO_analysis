#!/bin/bash
#SBATCH --job-name=recalc_Sj_coverage
#SBATCH --output=logs/recalc_Sj_coverage_%j.out
#SBATCH --error=logs/recalc_Sj_coverage_%j.err
#SBATCH --time=02:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=4

# Chargement des modules si nécessaire
module load samtools
module load minimap2

# Activation de l'environnement conda
source activate mahabio_env

# Chemins
CONTIGS="/shared/home/asandri/MAHABIO/data/Sj_contigs_more_than_300bp.fasta"
READS="/shared/home/asandri/MAHABIO/data/Sj.fastq.gz"
BAMOUT="/shared/home/asandri/MAHABIO_analysis/results/coverage/Sj_contigs_more_than_300bp_sorted.bam"
DEPTHOUT="/shared/home/asandri/MAHABIO_analysis/results/coverage/Sj_contigs_more_than_300bp_sorted.depth.txt"

echo "### Étape 1 : Alignement avec minimap2"
minimap2 -ax map-ont -t 4 "$CONTIGS" "$READS" | samtools view -@ 4 -bS - > temp_Sj.bam

echo "### Étape 2 : Tri du BAM"
samtools sort -@ 4 -o "$BAMOUT" temp_Sj.bam
rm temp_Sj.bam

echo "### Étape 3 : Indexation"
samtools index "$BAMOUT"

echo "### Étape 4 : Calcul de la couverture (depth)"
samtools depth "$BAMOUT" | awk '{sum[$1]+=$3} END {for (scaffold in sum) print scaffold, sum[scaffold]}' > "$DEPTHOUT"

echo "✅ Terminé. Fichiers générés :"
echo "  - $BAMOUT"
echo "  - ${BAMOUT}.bai"
echo "  - $DEPTHOUT"
