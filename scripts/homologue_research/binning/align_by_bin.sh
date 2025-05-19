#!/bin/bash
#SBATCH --job-name=align_bin_array
#SBATCH --output=logs/align_bin_%A_%a.out
#SBATCH --error=logs/align_bin_%A_%a.err
#SBATCH --time=00:30:00
#SBATCH --mem=3G
#SBATCH --cpus-per-task=2

source activate mahabio_env

# Dossiers
BASE="/shared/home/asandri/MAHABIO_analysis"
HITS_DIR="$BASE/results/binning/cmuA_hmmsearch/cmuA_hits_sequences"
ALIGN_DIR="$BASE/results/binning/cmuA_hmmsearch/alignment_phylo/per_bin"
REFS="/shared/home/asandri/MAHABIO/data/cmuaA_seq_prot_filtered.fasta"

mkdir -p "$ALIGN_DIR"

# Détection automatique des fichiers
FILES=($(ls $HITS_DIR/*_cmuA_hits.faa))
FILE=${FILES[$SLURM_ARRAY_TASK_ID]}
BASENAME=$(basename "$FILE" .faa)
CONCAT="$ALIGN_DIR/${BASENAME}_with_refs.faa"
ALIGNED="$ALIGN_DIR/${BASENAME}_aligned.faa"
TREE="$ALIGN_DIR/${BASENAME}_tree.nwk"

# Étapes
cat "$REFS" "$FILE" > "$CONCAT"
mafft --thread 2 --auto "$CONCAT" > "$ALIGNED"
FastTree -lg "$ALIGNED" > "$TREE"

echo "✅ Arbre généré pour $BASENAME"
