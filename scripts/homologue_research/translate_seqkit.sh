#!/bin/bash
#SBATCH --job-name=translate_seqkit
#SBATCH --output=logs/translate_seqkit_%j.out
#SBATCH --error=logs/seqkit_translate_seqkit_%j.err
#SBATCH --time=24:00:00
#SBATCH --mem=8G

source activate mahabio_env

mkdir -p /shared/home/asandri/MAHABIO/data/contigs_translated/contigs_translated_seqkit

for fasta in /shared/home/asandri/MAHABIO/data/*contigs*.fasta; do
    base=$(basename "$fasta" .fasta)
    seqkit translate --frame 6 -o "/shared/home/asandri/MAHABIO/data/contigs_translated/contigs_translated_seqkit/${base}_translated_raw.faa" "$fasta"
    #seqkit translate --frame 6 -p -i -o "/shared/home/asandri/MAHABIO/data/contigs_translated/contigs_translated_seqkit/${base}_translated.faa" "$fasta"

    i=1
    while read line; do
    if [[ $line == ">"* ]]; then
        echo "${line}_seq$i"
        ((i++))
    else
        echo "$line"
    fi
    done < "/shared/home/asandri/MAHABIO/data/contigs_translated/contigs_translated_seqkit/${base}_translated_raw.faa" > "/shared/home/asandri/MAHABIO/data/contigs_translated/contigs_translated_seqkit/${base}_translated.faa"

done

python header_sample_C_translated.py /shared/home/asandri/MAHABIO/data/contigs_translated/contigs_translated_seqkit/C_contigs_more_than_300bp_translated_raw.faa /shared/home/asandri/MAHABIO/data/contigs_translated/contigs_translated_seqkit/C_contigs_more_than_300bp_translated.faa

rm /shared/home/asandri/MAHABIO/data/contigs_translated/contigs_translated_seqkit/*_translated_raw.faa