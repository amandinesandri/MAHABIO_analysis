#!/usr/bin/env python3

import os
import pandas as pd

# Paramètres
hits_file = "/shared/home/asandri/MAHABIO_analysis/results/binning/Affiliation/16S_detection/16S_extracted/16S_hits.tsv"
output_file = os.path.join(os.path.dirname(hits_file), "shared_16S_contigs.tsv")

# Lecture
df = pd.read_csv(hits_file, sep='\t')

# Création d'un identifiant unique pour chaque bin
df['bin_id'] = df['binner'] + '/' + df['sample'] + '/' + df['bin_number'].astype(str)

# Compter le nombre de bins par contig
shared = (
    df.groupby('contig')['bin_id']
    .agg(lambda x: list(sorted(set(x))))
    .reset_index()
)

# Filtrer les contigs présents dans plus d’un bin
shared['n_bins'] = shared['bin_id'].apply(len)
shared = shared[shared['n_bins'] > 1]

# Sauvegarde
shared.to_csv(output_file, sep='\t', index=False)

print(f"✅ Contigs partagés détectés et enregistrés dans : {output_file}")
