#!/usr/bin/env python3

import os
import argparse
import pandas as pd

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--hits_file', required=True, help="Chemin vers le fichier 16S_hits.tsv")
    parser.add_argument('--sample', required=True, help="Nom de l'échantillon (e.g. C, Fk, Sj)")
    parser.add_argument('--output_file', required=True, help="Fichier TSV de sortie")
    args = parser.parse_args()

    df = pd.read_csv(args.hits_file, sep='\t')

    # Ne conserver que l’échantillon d'intérêt
    df = df[df['sample'] == args.sample]

    # Création d'un identifiant unique pour chaque bin
    df['bin_id'] = df['binner'] + '/' + df['sample'] + '/' + df['bin_number'].astype(str)

    # Grouper les contigs associés à plusieurs bins
    shared = (
        df.groupby('contig')['bin_id']
        .agg(lambda x: list(sorted(set(x))))
        .reset_index()
    )

    shared['n_bins'] = shared['bin_id'].apply(len)
    shared = shared[shared['n_bins'] > 1]

    shared.to_csv(args.output_file, sep='\t', index=False)
    print(f"✅ Contigs partagés détectés pour sample {args.sample} dans : {args.output_file}")

if __name__ == "__main__":
    main()
