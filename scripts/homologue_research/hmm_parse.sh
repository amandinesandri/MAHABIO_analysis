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

# # Extraction pour chaque .domtbl
# for domtbl in "$DOMTBL_DIR"/*hits.domtbl; do
#     base=$(basename "$domtbl" _hits.domtbl)
#     out_tsv="$OUT_DIR/${base}_summary.tsv"
#     echo "Parsing $domtbl → $out_tsv"
#     python "$SCRIPT" "$domtbl" "$out_tsv"
# done


BASE=/shared/home/asandri/MAHABIO_analysis
OUT_DIR=$BASE/results/hmm/summary
COMBINED=$OUT_DIR/combined_hits_summary.tsv

# Nettoyage initial
rm -f "$COMBINED"

# Récupérer la première ligne du premier fichier pour faire l'en-tête
head -n 1 "$OUT_DIR"/*_summary.tsv | head -n1 > "$COMBINED"

# Boucle sur tous les fichiers
for f in "$OUT_DIR"/*_summary.tsv; do
    tail -n +2 "$f" | awk -v sample_name="$(basename "$f" _summary.tsv)" 'BEGIN{OFS="\t"}{print sample_name, $0}' >> "$COMBINED"
done

echo "✅ Nouveau fichier combiné propre généré : $COMBINED"
