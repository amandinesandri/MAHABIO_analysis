import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from glob import glob

# Re-définir les répertoires après reset
identity_dir = "/mnt/data/identity_matrices"
heatmap_output_dir = "/mnt/data/heatmap_identity"
os.makedirs(heatmap_output_dir, exist_ok=True)

# Trouver tous les fichiers de matrices d'identité
identity_files = glob(os.path.join(identity_dir, "*_identity_matrix.tsv"))

generated_heatmaps = []

# Pour chaque fichier, tracer une heatmap
for identity_file in identity_files:
    try:
        df = pd.read_csv(identity_file, sep='\t', index_col=0)

        plt.figure(figsize=(10, 8))
        sns.heatmap(df, cmap="viridis", annot=False, fmt=".1f", square=True, 
                    cbar_kws={"label": "Identité (%)"}, linewidths=0.5)
        plt.title(os.path.basename(identity_file).replace("_identity_matrix.tsv", ""))
        plt.xticks(rotation=90, fontsize=6)
        plt.yticks(rotation=0, fontsize=6)
        plt.tight_layout()

        plot_path = os.path.join(
            heatmap_output_dir,
            os.path.basename(identity_file).replace("_identity_matrix.tsv", "_heatmap.png")
        )
        plt.savefig(plot_path, dpi=300)
        plt.close()
        generated_heatmaps.append(plot_path)
    except Exception as e:
        print(f"Erreur lors du traitement de {identity_file} : {e}")

