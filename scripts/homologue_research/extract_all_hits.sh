#!/bin/bash
#SBATCH --job-name=extract_all_hits
#SBATCH --output=logs/extract_all_hits_%j.out
#SBATCH --error=logs/extract_all_hits_%j.err
#SBATCH --time=00:20:00
#SBATCH --mem=4G

source activate mahabio_env

# Variables
BASE="/shared/home/asandri/MAHABIO_analysis"
SUMMARY_DIR="$BASE/results/hmm/summary"
FASTA_DIR="/shared/home/asandri/MAHABIO/data/contigs_translated/contigs_translated_seqkit"
EXTRACTED_DIR="$BASE/results/hmm/extracted_hits"
mkdir -p "$EXTRACTED_DIR"

# Boucle sur chaque fichier summary
for summary in "$SUMMARY_DIR"/*_summary.tsv; do
    sample=$(basename "$summary" _summary.tsv)
    faa="$FASTA_DIR/${sample}.faa"
    output="$EXTRACTED_DIR/${sample}_all_hits.faa"

    echo "Extraction de TOUS les hits pour $sample"
    cut -f1 "$summary" | tail -n +2 | seqkit grep -f - "$faa" > "$output"
done
