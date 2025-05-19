# Voici le petit script d'audit rapide demandé pour vérifier les GFF et comprendre
# pourquoi certains bins ont été comptés mais n'ont pas donné lieu à extraction.

import os
import pandas as pd

# === Paramètres ===
root_dir = "/shared/home/asandri/MAHABIO_analysis/results/binning/stats/16S_detection"
binners_samples = [
    ("semibin", "C_contigs_more_than_300bp"),
    ("vamb", "C_contigs_more_than_300bp"),
]

# === Audit ===
audit_records = []

for binner, sample in binners_samples:
    path = os.path.join(root_dir, binner, sample)
    if not os.path.exists(path):
        continue

    for file in os.listdir(path):
        if file.endswith(".gff"):
            file_path = os.path.join(path, file)
            size = os.path.getsize(file_path)

            n_16s = 0
            if size > 0:
                with open(file_path) as f:
                    for line in f:
                        if not line.startswith("#") and "16S" in line:
                            n_16s += 1

            conclusion = "✅ 16S détecté" if n_16s > 0 else "❌ Aucun 16S"
            audit_records.append({
                "binner": binner,
                "sample": sample,
                "gff_file": file,
                "gff_size (bytes)": size,
                "n_16S_annotations": n_16s,
                "conclusion": conclusion
            })

# === Résultats sous forme de tableau
audit_df = pd.DataFrame(audit_records)

import ace_tools as tools; tools.display_dataframe_to_user(name="Audit des fichiers GFF pour semibin_C et vamb_C", dataframe=audit_df)

audit_df.head()
