#!/bin/bash
#SBATCH --job-name=extract_matched
#SBATCH --output=logs/extract_matched_%j.out
#SBATCH --error=logs/extract_matched_%j.err
#SBATCH --time=01:00:00
#SBATCH --mem=4G

source activate mahabio_env

mkdir -p ../../results/hmm_matched

for domtbl in ../../results/hmm/*_hits.domtbl; do
    base=$(basename "$domtbl" _hits.domtbl)
    faa="/shared/home/asandri/MAHABIO/data/contigs_translated/contigs_translated_seqkit/${base}_filtered.faa"
    output="../../results/hmm_matched/${base}_matched.faa"

    python extract_matched_regions.py "$domtbl" "$faa" "$output"
done
