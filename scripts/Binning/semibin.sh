#!/bin/bash
#SBATCH --job-name=semibin
#SBATCH --output=logs/semibin_%A_%a.out
#SBATCH --error=logs/semibin_%A_%a.err
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=24:00:00

# Activation de l'environnement conda
source activate mahabio_env

# Définition des fichiers
contig_files=("/shared/home/asandri/MAHABIO/data/C_contigs_more_than_300bp.fasta" "/shared/home/asandri/MAHABIO/data/Sj_contigs_more_than_300bp.fasta" "/shared/home/asandri/MAHABIO/data/Fk_contigs_more_than_300bp.fasta")
bam_files=("/shared/home/asandri/MAHABIO_analysis/results/coverage/C_contigs_more_than_300bp.bam" "/shared/home/asandri/MAHABIO_analysis/results/coverage/Sj_contigs_more_than_300bp.bam" "/shared/home/asandri/MAHABIO_analysis/results/coverage/Fk_contigs_more_than_300bp.bam")

OUTDIR="/shared/home/asandri/MAHABIO_analysis/results/binning/semibin"
mkdir -p ${OUTDIR}

for i in "${!contig_files[@]}"; do
  CTG="${contig_files[$i]}"
  BAM="${bam_files[$i]}"
  BASENAME=$(basename "$CTG" .fasta)

  echo "Traitement de $BASENAME avec SemiBin..."

  semibin single_easy_bin \
    --contig "$CTG" \
    --bam "$BAM" \
    --outdir "${OUTDIR}/${BASENAME}" \
    --threads 16 \
    --environment soil
done
