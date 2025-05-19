#!/bin/bash
#SBATCH --job-name=cmuA_orf_positions
#SBATCH --output=logs/cmuA_orf_positions.out
#SBATCH --error=logs/cmuA_orf_positions.err
#SBATCH --time=00:05:00
#SBATCH --mem=1G

source activate mahabio_env

# mkdir -p /shared/home/asandri/MAHABIO_analysis/results/binning/cmuA_hmmsearch/visualization/

# awk -F'\t' 'BEGIN {OFS="\t"}
#     NR == 1 || ($1 == "semibin" && $2 == "Sj_contigs_more_than_300bp" && $3 == "20")
# ' /shared/home/asandri/MAHABIO_analysis/results/binning/cmuA_hmmsearch/cmuA_hits_summary_1e-5_50.tsv \
# > /shared/home/asandri/MAHABIO_analysis/results/binning/cmuA_hmmsearch/visualization/cmuA_bin20_hits.tsv

# # Exécution du script avec tes fichiers
# python extract_orf_positions.py \
#   /shared/home/asandri/MAHABIO_analysis/results/binning/cmuA_hmmsearch/visualization/cmuA_bin20_hits.tsv \
#   /shared/home/asandri/MAHABIO_analysis/results/binning/checkm/semibin/Sj_contigs_more_than_300bp/bins/SemiBin_20/genes.gff \
#   /shared/home/asandri/MAHABIO_analysis/results/binning/cmuA_hmmsearch/visualization/cmuA_positions_from_gff.tsv

python add_length_distance.py /shared/home/asandri/MAHABIO_analysis/results/binning/cmuA_hmmsearch/visualization/cmuA_positions_from_gff.tsv /shared/home/asandri/MAHABIO_analysis/results/binning/cmuA_hmmsearch/visualization/cmuA_positions_enriched.tsv
