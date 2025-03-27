#!/bin/bash
#SBATCH --job-name=hmm_parse
#SBATCH --output=logs/hmm_parse_%j.out
#SBATCH --error=logs/hmm_parse_%j.err
#SBATCH --time=01:00:00
#SBATCH --mem=8G

source activate mahabio_env

BASE=/shared/home/asandri/MAHABIO_analysis
DOMTBL_DIR=$BASE/results/hmm
OUT_DIR=$DOMTBL_DIR/summary
SCRIPT=$BASE/scripts/homologue_research/hmm_parse.py

mkdir -p "$OUT_DIR"

# Extraction pour chaque .domtbl
for domtbl in "$DOMTBL_DIR"/*hits.domtbl; do
    base=$(basename "$domtbl" _hits.domtbl)
    out_tsv="$OUT_DIR/${base}_summary.tsv"
    echo "Parsing $domtbl → $out_tsv"
    python "$SCRIPT" "$domtbl" "$out_tsv"
done

# Fusion de tous les fichiers
echo "Fusion des fichiers dans $OUT_DIR/combined_hits_summary.tsv"
head -n 1 "$OUT_DIR"/*_summary.tsv | head -n1 > "$OUT_DIR/combined_hits_summary.tsv"
for f in "$OUT_DIR"/*_summary.tsv; do
    tail -n +2 "$f" | awk -v sample="$(basename $f _summary.tsv)" '{print $0 "\t" sample}'
done >> "$OUT_DIR/combined_hits_summary.tsv"
