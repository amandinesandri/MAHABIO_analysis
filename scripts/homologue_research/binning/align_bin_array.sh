#!/bin/bash
#SBATCH --job-name=align_bin_array
#SBATCH --output=logs/align_bin_%A_%a.out
#SBATCH --error=logs/align_bin_%A_%a.err
#SBATCH --time=00:30:00
#SBATCH --mem=3G
#SBATCH --cpus-per-task=2

source activate mahabio_env
BASE="/shared/home/asandri/MAHABIO_analysis"

HITS_DIR="$BASE/results/binning/cmuA_hmmsearch/cmuA_hits_sequences"
N=$(ls $HITS_DIR/*_cmuA_hits.faa | wc -l)
sbatch --array=0-$(($N - 1)) $BASE/scripts/homologue_research/binning/align_by_bin.sh
