#!/bin/bash
#SBATCH --job-name=translate
#SBATCH --output=logs/translate_%j.out
#SBATCH --error=logs/translate_%j.err
#SBATCH --time=24:00:00
#SBATCH --mem=8G

module load conda
source deactivate
source activate mahabio_env

mkdir -p data/contigs_translated

for fasta in /shared/home/asandri/MAHABIO/data/*contigs*.fasta; do
    base=$(basename "$fasta" .fasta)
    transeq -sequence "$fasta" -outseq "data/contigs_translated/${base}_translated.faa" -frame 6
done