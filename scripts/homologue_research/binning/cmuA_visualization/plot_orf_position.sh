#!/bin/bash
#SBATCH --job-name=cmuA_plot
#SBATCH --output=logs/cmuA_plot.out
#SBATCH --error=logs/cmuA_plot.err
#SBATCH --time=00:10:00
#SBATCH --mem=2G

source activate mahabio_env
# Exécution du script Python
python plot_cmuA_orfs.py /shared/home/asandri/MAHABIO_analysis/results/binning/cmuA_hmmsearch/visualization/cmuA_positions_enriched.tsv
