#!/bin/bash
#SBATCH --job-name=cmuA_plot
#SBATCH --output=logs/cmuA_plot.out
#SBATCH --error=logs/cmuA_plot.err
#SBATCH --time=00:10:00
#SBATCH --mem=2G

source activate mahabio_env

mkdir -p /shared/home/asandri/MAHABIO_analysis/results/binning/cmuA_hmmsearch/visualization/

python plot_cmuA_position.py /shared/home/asandri/MAHABIO_analysis/results/binning/cmuA_hmmsearch/semibin/Sj_contigs_more_than_300bp/semibin_Sj_contigs_more_than_300bp_20_cmuA.tbl /shared/home/asandri/MAHABIO_analysis/results/binning/cmuA_hmmsearch/visualization/semibin_Sj_contigs_more_than_300bp_20_cmuA_hits_plot.png
