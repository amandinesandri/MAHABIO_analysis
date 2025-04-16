#!/bin/bash
#SBATCH --job-name=vamb
#SBATCH --output=logs/vamb_%A_%a.out
#SBATCH --error=logs/vamb_%A_%a.err
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --array=0-2

source activate vamb_env

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

i=$SLURM_ARRAY_TASK_ID

CTG="${contig_files[$i]}"
BAM="${bam_files[$i]}"
BASENAME=$(basename "$CTG" .fasta)
OUTDIR="/shared/home/asandri/MAHABIO_analysis/results/binning/vamb/${BASENAME}"
mkdir -p "$OUTDIR"

# Fichiers intermédiaires
TNF="${OUTDIR}/${BASENAME}.tnf"
RPKM="${OUTDIR}/${BASENAME}.rpkm"

echo "=== Génération des fichiers TNF et RPKM ==="
vamb --fasta "$CTG" --bamfiles "$BAM" --outdir "$OUTDIR" \
     --tnf "$TNF" --rpkm "$RPKM" --minfasta 200000 --threads 16

# Ou plus probablement :
# python -m vamb.tnf "$CTG" "$TNF"
# python -m vamb.bam2rpkm "$CTG" "$BAM" "$RPKM"

# Puis lancement du binning
echo "=== Lancement de VAMB classique ==="
vamb "$OUTDIR" "$TNF" "$RPKM" --minfasta 200000
