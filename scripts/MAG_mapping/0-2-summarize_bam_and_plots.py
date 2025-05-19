import os
import subprocess
import pandas as pd
import matplotlib.pyplot as plt

# -------------------
# PARAMÈTRES GÉNÉRAUX
# -------------------

bam_base = "/shared/home/asandri/MAHABIO_analysis/bam_files_full/"
output_dir = "/shared/home/asandri/MAHABIO_analysis/results/MAG_mappings/"
os.makedirs(output_dir, exist_ok=True)

methods = ["maxbin", "semibin", "vamb"]
samples = ["C", "Fk", "Sj"]
types = ["long", "short"]

SHORT_READS_COVERAGE_THRESHOLD = 30
LONG_READS_COVERAGE_THRESHOLD = 20

# -------------------
# EXTRACTION DES DONNÉES
# -------------------

data = []

print("▶️ Début de l'analyse des BAMs...")

for method in methods:
    for sample in samples:
        for readtype in types:
            bam_folder = f"{bam_base}/{method}/{sample}/{readtype}"
            if not os.path.exists(bam_folder):
                continue
            for bam_file in os.listdir(bam_folder):
                if not bam_file.endswith(".bam"):
                    continue
                bam_path = os.path.join(bam_folder, bam_file)
                mag_name = bam_file.replace("_long.bam", "").replace("_short.bam", "")

                flagstat = subprocess.check_output(f"samtools flagstat {bam_path}", shell=True).decode()
                mapped_reads = int(flagstat.split("\n")[4].split()[0])

                depth_output = subprocess.check_output(f"samtools depth -a {bam_path}", shell=True).decode().splitlines()
                depths = [int(line.split()[2]) for line in depth_output if len(line.split()) >= 3]
                mean_cov = sum(depths)/len(depths) if depths else 0

                data.append({
                    "Method": method,
                    "Sample": sample,
                    "MAG": mag_name,
                    "ReadType": readtype,
                    "MappedReads": mapped_reads,
                    "MeanCoverage": round(mean_cov, 2)
                })

print("✅ Extraction terminée.")

# -------------------
# CRÉATION DU TABLEAU COMPLET
# -------------------

df = pd.DataFrame(data)
print("✅ Aperçu rapide du DataFrame df :")
print(df.head())
pivot_df = df.pivot_table(index=["Method", "Sample", "MAG"], columns="ReadType", values=["MappedReads", "MeanCoverage"]).reset_index()
pivot_df.columns = ['Method', 'Sample', 'MAG',
                    'MappedReads_long', 'MappedReads_short',
                    'MeanCoverage_long', 'MeanCoverage_short']

# TRI décroissant par couverture short reads
pivot_df = pivot_df.sort_values(by="MeanCoverage_short", ascending=False)

pivot_df.to_csv(os.path.join(output_dir, "summary_MAGs_mapping.tsv"), sep="\t", index=False)
print(f"✅ Fichier complet trié sauvé : {output_dir}summary_MAGs_mapping.tsv")

# -------------------
# FILTRAGE DES BONS MAGs
# -------------------

def is_high_quality(row):
    return (row['MeanCoverage_short'] >= SHORT_READS_COVERAGE_THRESHOLD) and (row['MeanCoverage_long'] >= LONG_READS_COVERAGE_THRESHOLD)

pivot_df['HighQuality'] = pivot_df.apply(is_high_quality, axis=1)

high_quality_df = pivot_df[pivot_df['HighQuality'] == True]
high_quality_df.to_csv(os.path.join(output_dir, "high_quality_MAGs.tsv"), sep="\t", index=False)
print(f"✅ Fichier filtré sauvé : {output_dir}high_quality_MAGs.tsv")

# -------------------
# DOCUMENTATION DES SEUILS
# -------------------

with open(os.path.join(output_dir, "filtering_criteria.tsv"), "w") as f:
    f.write("Seuils de filtrage utilisés pour définir un MAG de haute qualité:\n")
    f.write(f"- Couverture moyenne short reads ≥ {SHORT_READS_COVERAGE_THRESHOLD}x\n")
    f.write(f"- Couverture moyenne long reads ≥ {LONG_READS_COVERAGE_THRESHOLD}x\n")
    f.write("\nJustification:\n")
    f.write("- ≥30x short reads garantit une détection fiable des gènes fonctionnels et une stabilité d'assemblage.\n")
    f.write("- ≥20x long reads compense les erreurs Nanopore et garantit une bonne continuité des contigs.\n")
print(f"✅ Critères documentés dans : {output_dir}filtering_criteria.tsv")

# -------------------
# PLOT DE QUALITÉ ROUGE/VERT
# -------------------

plt.figure(figsize=(20,10))
colors = pivot_df['HighQuality'].map({True: 'green', False: 'red'})

plt.bar(pivot_df['MAG'] + "_" + pivot_df['Sample'] + "_" + pivot_df['Method'],
        pivot_df['MeanCoverage_short'],
        color=colors,
        alpha=0.7)

plt.axhline(y=SHORT_READS_COVERAGE_THRESHOLD, color='blue', linestyle='--', label=f'Seuil {SHORT_READS_COVERAGE_THRESHOLD}x short reads')
plt.xticks(rotation=90, fontsize=5)
plt.ylabel("Mean Coverage Short Reads (X)")
plt.title("Qualité des MAGs basée sur la couverture short reads (Tri décroissant)")
plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(output_dir, "coverage_quality_plot.png"), dpi=300)
print(f"✅ Barplot qualité trié sauvé : {output_dir}coverage_quality_plot.png")

print("\n🏁 TOUT est terminé. Résultats dans", output_dir)