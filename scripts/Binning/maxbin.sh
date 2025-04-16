#!/bin/bash
#SBATCH --job-name=maxbin_binning
#SBATCH --output=logs/maxbin_%A_%a.out
#SBATCH --error=logs/maxbin_%A_%a.err
#SBATCH --time=08:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=8
#SBATCH --array=0-2

# 🧪 Active Conda dans le contexte SLURM (important !)
source activate mahabio_env

# 🔍 Vérification de l'environnement
echo "======== Environnement actif : $CONDA_DEFAULT_ENV ========"
echo "Bowtie2 path: $(which bowtie2)"
bowtie2 --version
echo "HMMER version:"
hmmsearch -h 2>&1 | grep HMMER
echo "FragGeneScan version:"
FragGeneScan -v

# 📁 Fichiers d'entrée (un par tâche)
contig_files=(
  "/shared/home/asandri/MAHABIO/data/C_contigs_more_than_300bp.fasta"
  "/shared/home/asandri/MAHABIO/data/Sj_contigs_more_than_300bp.fasta"
  "/shared/home/asandri/MAHABIO/data/Fk_contigs_more_than_300bp.fasta"
)

depth_files=(
  "/shared/home/asandri/MAHABIO_analysis/results/coverage/C_contigs_more_than_300bp_sorted.depth.txt"
  "/shared/home/asandri/MAHABIO_analysis/results/coverage/Sj_contigs_more_than_300bp_sorted.depth.txt"
  "/shared/home/asandri/MAHABIO_analysis/results/coverage/Fk_contigs_more_than_300bp_sorted.depth.txt"
)

reads=(
    "/shared/home/asandri/MAHABIO/data/C.fastq.gz"
    "/shared/home/asandri/MAHABIO/data/Sj.fastq.gz"
    "/shared/home/asandri/MAHABIO/data/Fk.fastq.gz"
)

# 🔢 Index de tâche courant
i=$SLURM_ARRAY_TASK_ID

CTG="${contig_files[$i]}"
DEPTH="${depth_files[$i]}"
READ="${reads[$i]}"
BASENAME=$(basename "$CTG" .fasta)

# 📂 Dossier de sortie
OUTDIR="/shared/home/asandri/MAHABIO_analysis/results/binning/maxbin"
SUBDIR="${OUTDIR}/${BASENAME}"
mkdir -p "$SUBDIR"

echo "======== Traitement de $BASENAME ========"

# # 🚀 Lancement de MaxBin2
# run_MaxBin.pl \
#   -contig "$CTG" \
#   -abund "$DEPTH" \
#   -out "${SUBDIR}/${BASENAME}_maxbin" \
#   -thread 8

# 🚀 Lancement de MaxBin2
run_MaxBin.pl \
  -contig "$CTG" \
  -reads "$READ" \
  -out "${SUBDIR}/${BASENAME}_maxbin" \
  -thread 8

echo "======== Fin du traitement de $BASENAME ========"
