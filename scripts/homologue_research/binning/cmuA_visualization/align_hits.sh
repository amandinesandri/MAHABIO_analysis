#!/bin/bash
#SBATCH --job-name=align_cmuA
#SBATCH --output=align_cmuA.out
#SBATCH --error=align_cmuA.err
#SBATCH --time=00:10:00
#SBATCH --mem=4G

source activate mahabio_env

muscle -in /shared/home/asandri/MAHABIO_analysis/results/binning/cmuA_hmmsearch/cmuA_hits_sequences/semibin_Sj_contigs_more_than_300bp_20_cmuA_hits.faa \
       -out /shared/home/asandri/MAHABIO_analysis/results/binning/cmuA_hmmsearch/visualization/cmuA_hits_aligned.faa
