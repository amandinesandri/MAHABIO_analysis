#!/bin/bash
#SBATCH --job-name=consensus
#SBATCH --output=logs/consensus_%j.out
#SBATCH --error=logs/consensus_%j.err
#SBATCH --time=00:30:00
#SBATCH --mem=4G

source activate mahabio_env

cons -sequence ../../results/hmm/tree/aligned_all_CmuA.fasta -outseq ../../results/hmm/cmuA_consensus.fasta -name cmuAcons