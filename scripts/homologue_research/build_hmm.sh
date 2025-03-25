#!/bin/bash
#SBATCH --job-name=build_hmm
#SBATCH --output=logs/build_hmm.out
#SBATCH --error=logs/build_hmm.err
#SBATCH --time=01:00:00
#SBATCH --mem=4G

module load conda
source activate mahabio_env

cd /shared/home/asandri/MAHABIO_analysis

mafft --auto /shared/home/asandri/MAHABIO/data/cmuaA_seq_prot.fasta > results/hmm/cmuA_aligned.fasta

hmmbuild results/hmm/cmuA.hmm results/hmm/cmuA_aligned.fasta

source deactivate