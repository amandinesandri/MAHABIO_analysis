import os
from Bio import Phylo
import matplotlib.pyplot as plt

# Dossier contenant les arbres
tree_dir = "/shared/home/asandri/MAHABIO_analysis/results/binning/cmuA_hmmsearch/alignment_phylo/per_bin"  # ← À ADAPTER
output_dir = os.path.join(tree_dir, "png")
os.makedirs(output_dir, exist_ok=True)

tree_files = [f for f in os.listdir(tree_dir) if f.endswith("_tree.nwk")]

for tree_file in tree_files:
    tree_path = os.path.join(tree_dir, tree_file)
    tree = Phylo.read(tree_path, "newick")

    # Création de la figure
    fig = plt.figure(figsize=(10, 12))
    ax = fig.add_subplot(1, 1, 1)
    Phylo.draw(tree, do_show=False, axes=ax)
    ax.set_title(tree_file.replace("_tree.nwk", ""), fontsize=10)

    # Sauvegarde
    output_path = os.path.join(output_dir, tree_file.replace(".nwk", ".png"))
    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    plt.close(fig)

    print(f"✅ Exporté : {output_path}")
