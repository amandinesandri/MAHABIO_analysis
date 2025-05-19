import pandas as pd

# Charger les données brutes
df = pd.read_csv("/shared/home/asandri/MAHABIO_analysis/results/binning/cmuA_hmmsearch/alignment_phylo/per_bin/cmuA_best_hits_summary.tsv", sep="\t")

# Vérification qu'on a bien plusieurs bin_seq/closest_ref par bin
# Puis extraire pour chaque bin la ligne avec la plus petite distance
grouped = []
for bin_name, group in df.groupby("bin"):
    best_row = group.loc[group["distance"].idxmin()]
    grouped.append(best_row)

best_hits_df = pd.DataFrame(grouped)

# Sauvegarde du fichier corrigé
output_path = "/shared/home/asandri/MAHABIO_analysis/results/binning/cmuA_hmmsearch/alignment_phylo/per_bin/cmuA_best_hit_per_bin.tsv"
best_hits_df.to_csv(output_path, sep="\t", index=False)
output_path
