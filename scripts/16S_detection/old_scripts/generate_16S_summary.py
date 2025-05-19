#!/usr/bin/env python3

import os
import pandas as pd

# Paramètres
hits_file = "/shared/home/asandri/MAHABIO_analysis/results/binning/Affiliation/16S_detection/16S_extracted/16S_hits.tsv"
output_file = os.path.join(os.path.dirname(hits_file), "16S_hits_summary.tsv")

# Lecture du fichier
df = pd.read_csv(hits_file, sep='\t')

# Création du tableau résumé
summary = (
    df.groupby(['binner', 'sample', 'bin_number'])
    .agg(
        n_hits_16S=('Name', 'count'),
        contigs=('contig', lambda x: ', '.join(sorted(set(x))))
    )
    .reset_index()
)

# Sauvegarde
summary.to_csv(output_file, sep='\t', index=False)

print(f"✅ Résumé généré : {output_file}")
