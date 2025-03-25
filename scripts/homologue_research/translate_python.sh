#!/bin/bash
#SBATCH --job-name=translate_python
#SBATCH --output=logs/translate_python_%j.out
#SBATCH --error=logs/seqkit_translate_python_%j.err
#SBATCH --time=24:00:00
#SBATCH --mem=8G

source activate mahabio_env

mkdir -p /shared/home/asandri/MAHABIO/data/contigs_translated/contigs_translated_seqkit/contigs_translated_python

for fasta in /shared/home/asandri/MAHABIO/data/*contigs*.fasta; do
    base=$(basename "$fasta" .fasta)
    python /shared/home/asandri/MAHABIO_analysis/scripts/homologue_research/translate_6frames.py "$fasta" "/shared/home/asandri/MAHABIO/data/contigs_translated/contigs_translated_python/${base}_translated.faa"
done