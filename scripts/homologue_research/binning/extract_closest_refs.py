import os
from ete3 import Tree
import pandas as pd

# Configuration
TREE_DIR = "/shared/home/asandri/MAHABIO_analysis/results/binning/cmuA_hmmsearch/alignment_phylo/per_bin"
REF_KEYWORDS = ["MB2", "PA1", "BJ001", "ref_", "LT674", "AY439", "AY934", "GAF", "GAG", "GAI"]  # adapte si besoin
OUTFILE = "/shared/home/asandri/MAHABIO_analysis/results/binning/cmuA_hmmsearch/alignment_phylo/per_bin/cmuA_best_hits_summary.tsv"

summary = []

def is_reference(name):
    return any(key in name for key in REF_KEYWORDS)

# Analyse de tous les arbres
for file in os.listdir(TREE_DIR):
    if file.endswith("_tree.nwk"):
        path = os.path.join(TREE_DIR, file)
        tree = Tree(path)

        bin_seqs = [leaf.name for leaf in tree if not is_reference(leaf.name)]
        ref_seqs = [leaf.name for leaf in tree if is_reference(leaf.name)]

        if not bin_seqs or not ref_seqs:
            continue

        for bseq in bin_seqs:
            node_b = tree.search_nodes(name=bseq)[0]
            closest_ref = None
            min_dist = float("inf")

            for rseq in ref_seqs:
                node_r = tree.search_nodes(name=rseq)[0]
                dist = tree.get_distance(node_b, node_r)
                if dist < min_dist:
                    min_dist = dist
                    closest_ref = rseq

            summary.append({
                "bin": file.replace("_tree.nwk", ""),
                "bin_seq": bseq,
                "closest_ref": closest_ref,
                "distance": round(min_dist, 4)
            })

# Export
df = pd.DataFrame(summary)
df.to_csv(OUTFILE, sep="\t", index=False)
print(f"✅ Résultats enregistrés dans {OUTFILE}")
