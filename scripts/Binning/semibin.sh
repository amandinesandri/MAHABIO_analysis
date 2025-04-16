#!/bin/bash
#SBATCH --job-name=semibin
#SBATCH --output=logs/semibin_%A_%a.out
#SBATCH --error=logs/semibin_%A_%a.err
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --array=0-2

# Activation de l'environnement conda
source activate mahabio_env

# Paramètres
contig_files=(
  "/shared/home/asandri/MAHABIO/data/C_contigs_more_than_300bp.fasta"
  "/shared/home/asandri/MAHABIO/data/Sj_contigs_more_than_300bp.fasta"
  "/shared/home/asandri/MAHABIO/data/Fk_contigs_more_than_300bp.fasta"
)

bam_files=(
  "/shared/home/asandri/MAHABIO_analysis/results/coverage/C_contigs_more_than_300bp_sorted.bam"
  "/shared/home/asandri/MAHABIO_analysis/results/coverage/Sj_contigs_more_than_300bp_sorted.bam"
  "/shared/home/asandri/MAHABIO_analysis/results/coverage/Fk_contigs_more_than_300bp_sorted.bam"
)

# Index de tâche SLURM
i=$SLURM_ARRAY_TASK_ID

CTG="${contig_files[$i]}"
BAM="${bam_files[$i]}"
BASENAME=$(basename "$CTG" .fasta)
OUTDIR="/shared/home/asandri/MAHABIO_analysis/results/binning/semibin/${BASENAME}"

mkdir -p "$OUTDIR"

echo "=== Traitement de $BASENAME avec SemiBin ==="

semibin single_easy_bin \
  --contig "$CTG" \
  --bam "$BAM" \
  --outdir "$OUTDIR" \
  --threads 16 \
  --environment soil

echo "=== Fin du traitement de $BASENAME ==="
