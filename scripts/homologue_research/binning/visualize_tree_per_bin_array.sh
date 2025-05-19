#!/bin/bash
#SBATCH --job-name=visualize_tree
#SBATCH --output=logs/visualize_tree%A_%a.out
#SBATCH --error=logs/visualize_tree_%A_%a.err
#SBATCH --time=00:30:00
#SBATCH --mem=3G
#SBATCH --cpus-per-task=2

source activate mahabio_env

python /shared/home/asandri/MAHABIO_analysis/scripts/homologue_research/binning/visualize_tree_per_bin.py