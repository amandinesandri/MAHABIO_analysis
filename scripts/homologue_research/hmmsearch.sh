#!/bin/bash
#SBATCH --job-name=hmmsearch
#SBATCH --output=logs/hmmsearch_%j.out
#SBATCH --error=logs/hmmsearch_%j.err
#SBATCH --time=02:00:00
#SBATCH --mem=16G


module load hmmer


for protfile in /shared/home/asandri/MAHABIO/data/contigs_translated/contigs_translated_seqkit/*.faa; do
    base=$(basename "$protfile" .faa)
    filtered="/shared/home/asandri/MAHABIO/data/contigs_translated/contigs_translated_seqkit/${base}_filtered.faa"

    # Filtrage longueur max
    seqkit seq -m 30 -M 100000 "$protfile" > "$filtered"

    # hmmsearch sur le fichier filtré
    hmmsearch --tblout ../../results/hmm/"${base}_filtered_hits.tbl" --domtblout ../../results/hmm/"${base}_filtered_hits.domtbl ../../results/hmm/cmuA.hmm "$filtered"

done