#!/bin/bash
#SBATCH --job-name=Alignment
#SBATCH --output=logs/Alignment_%j.out
#SBATCH --error=logs/Alignment_%j.err
#SBATCH --time=00:10:00
#SBATCH --mem=2G
#SBATCH --cpus-per-task=1

source activate mahabio_env

# Variables
BASE="/shared/home/asandri/MAHABIO_analysis"
HITS_DIR="$BASE/results/binning/cmuA_hmmsearch/cmuA_hits_sequences"
ALIGN_DIR="$BASE/results/binning/cmuA_hmmsearch/alignment_phylo"
mkdir -p "$ALIGN_DIR"

# # 1. Fusionner toutes les séquences hits
# cat "$HITS_DIR"/*.faa > "$ALIGN_DIR/all_cmuA_hits_raw.faa"
# echo "✅ Fusion brute terminée."

# # 2. Corriger les headers pour qu'ils soient absolument uniques
# awk '/^>/{print ">seq_" NR; next} {print}' "$ALIGN_DIR/all_cmuA_hits_raw.faa" > "$ALIGN_DIR/all_cmuA_hits_unique.faa"
# echo "✅ Headers renommés en seq_1, seq_2, etc."

# 3. MAFFT Alignement
mafft --thread 4 --auto "$ALIGN_DIR/all_cmuA_hits_unique.faa" > "$ALIGN_DIR/all_cmuA_hits_aligned.faa"
echo "✅ Alignement MAFFT terminé."

# 4. FastTree arbre phylogénétique
FastTree -lg "$ALIGN_DIR/all_cmuA_hits_aligned.faa" > "$ALIGN_DIR/all_cmuA_hits_tree.nwk"
echo "✅ Arbre FastTree terminé."
