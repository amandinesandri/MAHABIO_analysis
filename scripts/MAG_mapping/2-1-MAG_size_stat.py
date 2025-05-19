import pandas as pd

# Fichier d'entrée
input_file = "/shared/home/asandri/MAHABIO_analysis/results/MAG_mapping/MAG_sizes.tsv"

# Dossier de sortie
output_dir = "/shared/home/asandri/MAHABIO_analysis/results/MAG_mapping/"

# Charger les tailles des MAGs
df = pd.read_csv(input_file, sep="\t")

# Ajouter la classification par taille
df["Size_Class"] = df["Size_bp"].apply(lambda x: ">500kb" if x >= 500000 else "<500kb")

# Filtrer les MAGs >500kb
filtered_df = df[df["Size_bp"] >= 500000]
filtered_df.to_csv(f"{output_dir}/Bins_filtered_over_500kb.tsv", sep="\t", index=False)

# Calcul des statistiques par couple Sample/binner
grouped = df.groupby(["Sample", "binner"])

stats_list = []
for (sample, binner), group in grouped:
    total = len(group)
    mean_total = group["Size_bp"].mean()
    high = group[group["Size_bp"] >= 500000]
    low = group[group["Size_bp"] < 500000]
    stats_list.append({
        "Sample": sample,
        "Binner": binner,
        "Total_Bins": total,
        "Total_Bin_Size_Mean(bp)": round(mean_total, 2),
        "Bins_>500kb_Count": len(high),
        "Bins_>500kb_MeanSize(bp)": round(high["Size_bp"].mean(), 2) if not high.empty else 0,
        "Bins_<500kb_Count": len(low),
        "Bins_<500kb_MeanSize(bp)": round(low["Size_bp"].mean(), 2) if not low.empty else 0
    })

# Création du DataFrame résumé
stats_df = pd.DataFrame(stats_list)
stats_df.to_csv(f"{output_dir}/MAGs_size_summary.tsv", sep="\t", index=False)

print("✅ Statistiques résumées générées dans MAGs_size_summary.tsv")
print("✅ Fichier filtré généré dans Bins_filtered_>500kb.tsv")

# Représentation graphique
import seaborn as sns
import matplotlib.pyplot as plt

plt.figure(figsize=(12, 6))
sns.histplot(data=df, x="Size_Mbp", hue="binner", element="step", bins=50, stat="count", multiple="stack", palette="Set2")
plt.title("Distribution of MAG Sizes (in Mbp) by Binning Method")
plt.xlabel("MAG Size (Mbp)")
plt.ylabel("Number of MAGs")
plt.grid(axis="y")
plt.tight_layout()
plt.savefig(f"{output_dir}/Distribution_taille_MAGs.png", dpi=300)
plt.close()
