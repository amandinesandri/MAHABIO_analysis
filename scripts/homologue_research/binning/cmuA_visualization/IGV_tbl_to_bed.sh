#!/bin/bash
#SBATCH --job-name=tbl2bed
#SBATCH --output=tbl2bed.out
#SBATCH --error=tbl2bed.err
#SBATCH --time=00:05:00
#SBATCH --mem=1G

source activate mahabio_env

python IGV_tbl_to_bed.py /shared/home/asandri/MAHABIO_analysis/results/binning/cmuA_hmmsearch/semibin/Sj_contigs_more_than_300bp/semibin_Sj_contigs_more_than_300bp_20_cmuA.tbl /shared/home/asandri/MAHABIO_analysis/results/binning/cmuA_hmmsearch/visualization/cmuA_hits.bed
