import pandas as pd
import os
from glob import glob

# === Paramètres ===
base_path = "/shared/home/asandri/MAHABIO_analysis/results/binning/"
binners = ["maxbin", "semibin", "vamb"]
samples = ["C_contigs_more_than_300bp", "Fk_contigs_more_than_300bp", "Sj_contigs_more_than_300bp"]
output_file = "/shared/home/asandri/MAHABIO_analysis/results/MAG_mapping/merged_MAG_info.tsv"

# === Initialiser une liste pour tout stocker ===
merged_data = []

# === Boucle sur binner/sample ===
for binner in binners:
    for sample in samples:
        # chemins dRep et GTDBtk
        drep_info_path = os.path.join(base_path, "dRep", binner, sample, "data_tables", "genomeInformation.csv")
        gtdbtk_summary_path = os.path.join(base_path, "GTDBtk", binner, sample, "gtdbtk.bac120.summary.tsv")
        
        # Vérifier que les fichiers existent
        if not os.path.isfile(drep_info_path) or not os.path.isfile(gtdbtk_summary_path):
            continue
        
        # Charger les fichiers
        drep_df = pd.read_csv(drep_info_path)
        gtdbtk_df = pd.read_csv(gtdbtk_summary_path, sep='\t')
        
        # Nettoyer colonnes intéressantes
        drep_df = drep_df[['genome', 'completeness', 'contamination', 'Genome size (bp)', 'N50 (contigs)', 'contig_num']]
        gtdbtk_df = gtdbtk_df[['user_genome', 'classification']]
        
        # Merge sur les IDs des MAGs
        merged = pd.merge(drep_df, gtdbtk_df, how='left', left_on='genome', right_on='user_genome')
        
        # Extraire Famille et Genre
        families = []
        genres = []
        for tax in merged['classification'].fillna('d__Bacteria;p__;c__;o__;f__;g__').tolist():
            parts = tax.split(';')
            family = parts[4][3:] if len(parts) > 4 else 'Unknown'
            genus = parts[5][3:] if len(parts) > 5 else 'Unknown'
            families.append(family)
            genres.append(genus)
        
        merged['Family'] = families
        merged['Genus'] = genres
        
        # Ajouter infos supplémentaires
        merged['Binner'] = binner
        merged['Sample'] = sample
        merged['Source'] = f"{binner}/{sample}"
        
        # Garder colonnes importantes
        merged = merged[['genome', 'Binner', 'Sample', 'Source', 'completeness', 'contamination', 
                         'Genome size (bp)', 'N50 (contigs)', 'contig_num', 'Family', 'Genus']]
        
        merged_data.append(merged)

# === Concaténer tous les MAGs ensemble ===
final_df = pd.concat(merged_data, ignore_index=True)

# === Sauvegarder en TSV ===
final_df.to_csv(output_file, sep='\t', index=False)

print(f"Fichier fusionné créé : {output_file}")
