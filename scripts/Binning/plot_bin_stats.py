
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# === Charger les données
df = pd.read_csv("/shared/home/asandri/MAHABIO_analysis/results/binning/stats/all_binners_stats.tsv", sep="\t")

# === Nettoyage des données
df["Completeness"] = pd.to_numeric(df["Completeness"], errors="coerce")
df["Contamination"] = pd.to_numeric(df["Contamination"], errors="coerce")

# === Plot 1 : boxplot de la complétude
plt.figure(figsize=(8, 5))
sns.boxplot(x="binner", y="Completeness", data=df)
plt.title("completude Distribution")
plt.ylabel("Completude (%)")
plt.xlabel("Binner")
plt.tight_layout()
plt.savefig("/shared/home/asandri/MAHABIO_analysis/results/binning/stats/boxplot_completeness.png")
plt.close()

# === Plot 2 : boxplot de la contamination
plt.figure(figsize=(8, 5))
sns.boxplot(x="binner", y="Contamination", data=df)
plt.title("contamination distribution")
plt.ylabel("Contamination (%)")
plt.xlabel("Binner")
plt.tight_layout()
plt.savefig("/shared/home/asandri/MAHABIO_analysis/results/binning/stats/boxplot_contamination.png")
plt.close()

# === Plot 3 : nombre de bins récupérés
bin_counts = df["binner"].value_counts().reset_index()
bin_counts.columns = ["binner", "n_bins"]

plt.figure(figsize=(6, 4))
sns.barplot(x="binner", y="n_bins", data=bin_counts)
plt.title("Total Bins found")
plt.ylabel("Bin counts")
plt.xlabel("Binner")
plt.tight_layout()
plt.savefig("/shared/home/asandri/MAHABIO_analysis/results/binning/stats/bins_per_binner.png")
plt.close()

print("✅ Plots enregistrés dans results/stats/")