#!/bin/bash
#SBATCH --job-name=coverage
#SBATCH --output=logs/coverage_%j.out
#SBATCH --error=logs/coverage_%j.err
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
# CTG="/shared/home/asandri/MAHABIO/data/C_contigs_more_than_300bp.fasta"
# R1="/shared/home/asandri/MAHABIO/data/C_R1.fastq.gz"
# R2="/shared/home/asandri/MAHABIO/data/C_R2.fastq.gz"

# BASENAME=$(basename "$CTG" .fasta)

# echo "### Indexation de $CTG avec BWA"
# bwa index "$CTG"

# echo "### Alignement des reads de $BASENAME"
# bwa mem -t 8 "$CTG" "$R1" "$R2" | samtools view -bS - > "${BASENAME}.bam"

# echo "### Tri et indexation"
# samtools sort "${BASENAME}.bam" -o "${BASENAME}_sorted.bam"
# samtools index "${BASENAME}_sorted.bam"

# # Nettoyage si besoin
# #mv "${BASENAME}_sorted.bam" "${BASENAME}.bam"
# #mv "${BASENAME}_sorted.bam.bai" "${BASENAME}.bam.bai"
# #rm "${BASENAME}.bam"




# contigs=("/shared/home/asandri/MAHABIO/data/Sj_contigs_more_than_300bp.fasta" "/shared/home/asandri/MAHABIO/data/Fk_contigs_more_than_300bp.fasta")
# reads=("/shared/home/asandri/MAHABIO/data/Sj.fastq.gz" "/shared/home/asandri/MAHABIO/data/Fk.fastq.gz")

# for i in "${!contigs[@]}"; do
#     CTG="${contigs[$i]}"
#     READS="${reads[$i]}"
#     BASENAME=$(basename "$CTG" .fasta)

#     echo "### Alignement de $READS sur $CTG avec minimap2"
#     minimap2 -ax map-ont -t 8 "$CTG" "$READS" | samtools view -bS - > "${BASENAME}.bam"

#     echo "### Tri et indexation"
#     samtools sort "${BASENAME}.bam" -o "${BASENAME}_sorted.bam"
#     samtools index "${BASENAME}_sorted.bam"

#     #mv "${BASENAME}_sorted.bam" "${BASENAME}.bam"
#     #mv "${BASENAME}_sorted.bam.bai" "${BASENAME}.bam.bai"
#     #rm "${BASENAME}.bam"
# done

bam_files=("/shared/home/asandri/MAHABIO_analysis/results/coverage/Sj_contigs_more_than_300bp_sorted.bam" "/shared/home/asandri/MAHABIO_analysis/results/coverage/Fk_contigs_more_than_300bp_sorted.bam" "/shared/home/asandri/MAHABIO_analysis/results/coverage/C_contigs_more_than_300bp_sorted.bam")


for BAM in "${bam_files[@]}"; do


    BASENAME=$(basename "$BAM" .bam)
    echo "Traitement de $BAM..."
    
    # Indexation si besoin
    if [ ! -f "${BAM}.bai" ]; then
        samtools index "$BAM"
    fi

    # Génération du fichier depth
    samtools depth "$BAM" | awk '{sum[$1]+=$3} END {for (scaffold in sum) print scaffold, sum[scaffold]}' > "${BASENAME}.depth.txt"
done