import os
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

# Dossier contenant les matrices d'identité (au format TSV)
identity_dir = "/shared/home/asandri/MAHABIO_analysis/results/binning/Affiliation/16S_evaluation/identity"

# Liste des fichiers TSV
tsv_files = [f for f in os.listdir(identity_dir) if f.endswith("_identity_matrix.tsv")]

# Création d'un sous-dossier pour les heatmaps
heatmap_dir = os.path.join(identity_dir, "heatmaps")
os.makedirs(heatmap_dir, exist_ok=True)

# Chargement et visualisation de chaque matrice
for file in tsv_files:
    path = os.path.join(identity_dir, file)
    df = pd.read_csv(path, sep="\t", index_col=0)
    
    plt.figure(figsize=(8, 6))
    sns.heatmap(df, annot=True, fmt=".1f", cmap="viridis", cbar_kws={'label': 'Identité (%)'})
    plt.title(file.replace("_identity_matrix.tsv", "").replace("_", " "))
    plt.xticks(rotation=45, ha="right")
    plt.yticks(rotation=0)
    plt.tight_layout()
    
    plot_path = os.path.join(heatmap_dir, file.replace(".tsv", ".png"))
    plt.savefig(plot_path, dpi=300)
    plt.close()

# Liste des heatmaps générées
os.listdir(heatmap_dir)
