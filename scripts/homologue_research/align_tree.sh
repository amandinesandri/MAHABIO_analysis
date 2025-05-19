#!/bin/bash
#SBATCH --job-name=align_tree_all_hits
#SBATCH --output=logs/align_tree_all_hits_%A_%a.out
#SBATCH --error=logs/align_tree_all_hits_%A_%a.err
#SBATCH --time=01:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=4
#SBATCH --array=0-2    # <-- à adapter au nombre d'échantillons

source activate mahabio_env

# === Variables générales
BASE="/shared/home/asandri/MAHABIO_analysis"
EXTRACTED_DIR="$BASE/results/hmm/extracted_hits"
ALIGN_DIR="$BASE/results/hmm/alignment_phylo"
mkdir -p "$ALIGN_DIR"

# === Lister les fichiers
FILES=($(ls "$EXTRACTED_DIR"/*_all_hits.faa))

# === Choisir le fichier correspondant à l'index SLURM_ARRAY_TASK_ID
faa="${FILES[$SLURM_ARRAY_TASK_ID]}"
sample=$(basename "$faa" _all_hits.faa)

aligned="$ALIGN_DIR/${sample}_aligned.faa"
tree="$ALIGN_DIR/${sample}_tree.nwk"

echo "📈 Alignement avec MAFFT pour $sample..."
mafft --thread 4 --auto "$faa" > "$aligned"

echo "🌳 Construction de l'arbre phylogénétique avec FastTree pour $sample..."
FastTree -lg "$aligned" > "$tree"

echo "✅ $sample : alignement et arbre phylo terminés."
