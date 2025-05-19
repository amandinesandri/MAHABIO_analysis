# Chemin d'entrée du fichier de tailles
input_file = "/shared/home/asandri/MAHABIO_analysis/results/MAG_mappings/MAG_sizes.tsv"

# Dossier de sortie
output_dir = "/shared/home/asandri/MAHABIO_analysis/results/MAG_mappings/"

# Charger le fichier
df = pd.read_csv(input_file, sep="\t")

# Séparer selon la taille 500000 bp (0.5 Mbp)
high_size = df[df["Size_bp"] >= 500000]
low_size = df[df["Size_bp"] < 500000]

# Sauvegarder
high_size.to_csv(f"{output_dir}/high_size_MAGs.tsv", sep="\t", index=False)
low_size.to_csv(f"{output_dir}/low_size_MAGs.tsv", sep="\t", index=False)

print(f"✅ {len(high_size)} MAGs >500kb sauvegardés dans high_size_MAGs.tsv")
print(f"✅ {len(low_size)} MAGs <500kb sauvegardés dans low_size_MAGs.tsv")