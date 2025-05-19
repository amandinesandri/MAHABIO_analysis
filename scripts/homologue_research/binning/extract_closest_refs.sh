#!/bin/bash
#SBATCH --job-name=extract_closest_refs
#SBATCH --output=logs/extract_closest_refs_%A_%a.out
#SBATCH --error=logs/extract_closest_refs_%A_%a.err
#SBATCH --time=00:30:00
#SBATCH --mem=3G
#SBATCH --cpus-per-task=2

source activate mahabio_env

#python extract_closest_refs.py
python /shared/home/asandri/MAHABIO_analysis/scripts/homologue_research/binning/cmuA_best_hits_per_bin.py