#!/bin/bash
#SBATCH --job-name=reroot_tree
#SBATCH --output=logs/reroot_tree_%j.out
#SBATCH --error=logs/reroot_tree_%j.err
#SBATCH --time=00:20:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1

source activate mahabio_env

# Variables
BASE="/shared/home/asandri/MAHABIO_analysis"
TREE_IN="$BASE/results/binning/cmuA_hmmsearch/alignment_phylo/all_cmuA_hits_tree.nwk"
TREE_OUT="$BASE/results/binning/cmuA_hmmsearch/alignment_phylo/all_cmuA_hits_tree_rerooted.nwk"

# Choisis ici la séquence que tu veux utiliser comme racine
ROOT_NAME="seq_1"  # à adapter selon ton mapping si besoin !

# Utilisation de ETE3 pour rerooter
python - << EOF
from ete3 import Tree

# Charger l'arbre
t = Tree("$TREE_IN", format=1)

# Chercher le noeud racine
root_node = t.search_nodes(name="$ROOT_NAME")[0]

# Re-rooter
t.set_outgroup(root_node)

# Sauvegarder
t.write(format=1, outfile="$TREE_OUT")
print("✅ Arbre rerooté sauvegardé ici : $TREE_OUT")
EOF
