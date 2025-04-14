#!/bin/bash
#SBATCH --job-name=maxbin_binning
#SBATCH --output=logs/maxbin_%A_%a.out
#SBATCH --error=logs/maxbin_%A_%a.err
#SBATCH --time=08:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=8
#SBATCH --array=0-2

source activate mahabio_env


# Liste des fichiers d'entrée (fichier contigs et fichier couverture)
contig_files=("/shared/home/asandri/MAHABIO/data/C_contigs_more_than_300bp.fasta" "/shared/home/asandri/MAHABIO/data/Sj_contigs_more_than_300bp.fasta" "/shared/home/asandri/MAHABIO/data/Fk_contigs_more_than_300bp.fasta")
depth_files=("/shared/home/asandri/MAHABIO_analysis/results/coverage/C_contigs_more_than_300bp.depth.txt" "/shared/home/asandri/MAHABIO_analysis/results/coverage/Sj_contigs_more_than_300bp.depth.txt" "/shared/home/asandri/MAHABIO_analysis/results/coverage/Fk_contigs_more_than_300bp.depth.txt")


# Création du dossier de sortie

OUTDIR="/shared/home/asandri/MAHABIO_analysis/results/binning/maxbin"

mkdir -p ${OUTDIR}
cd ${OUTDIR}

for i in "${!contig_files[@]}"; do
  CTG="${contig_files[$i]}"
  DEPTH="${depth_files[$i]}"
  BASENAME=$(basename "$CTG" .fasta)  # Lancement de MaxBin2
  run_MaxBin.pl \
    -contig "$CTG" \
    -abund "$DEPTH" \
    -out "${OUTDIR}/${BASENAME}/${BASENAME}_maxbin" \
    -thread 8
done
