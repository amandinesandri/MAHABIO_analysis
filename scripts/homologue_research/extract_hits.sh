#!/bin/bash
#SBATCH --job-name=extract_hits
#SBATCH --output=logs/extract_hits_%j.out
#SBATCH --error=logs/extract_hits_%j.err
#SBATCH --time=01:00:00
#SBATCH --mem=4G

source activate mahabio_env

mkdir -p ../../results/hmm/hits

for tbl in ../../results/hmm/*_filtered_hits.tbl; do
    base=$(basename "$tbl" _filtered_hits.tbl)
    faa="/shared/home/asandri/MAHABIO/data/contigs_translated/contigs_translated_seqkit/${base}_filtered.faa"
    output="../../results/hmm/hits/${base}_hits.faa"

    python extract_hmm_hits.py "$tbl" "$faa" "$output"
done