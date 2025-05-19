#!/usr/bin/env python

# === SCRIPT BARRNAP SANS RISQUE DE COLLISION ET CORRECTION DU SUMMARY FINAL ===

import os
import subprocess
import pandas as pd
import argparse
import gzip
import shutil

# === Parser ===
parser = argparse.ArgumentParser()
parser.add_argument('--binning_root', required=True, help="Racine du dossier de binning")
parser.add_argument('--output_root', required=True, help="Racine du dossier pour stocker FASTA+GFF et summary")
args = parser.parse_args()

binners_info = {
    'semibin': {'pattern': 'output_bins', 'ext': '.fa.gz'},
    'maxbin': {'pattern': '', 'ext': '.fasta'},
    'vamb': {'pattern': 'bins', 'ext': '.fna'},
}

# === Fonction pour lancer barrnap ===
def run_barrnap_on_file(fasta_in, gff_out):
    cmd = f"barrnap --kingdom bac {fasta_in} > {gff_out}"
    result = subprocess.run(cmd, shell=True)
    if result.returncode != 0:
        print(f"⚠️ Barrnap a échoué pour {fasta_in}")
        return False
    return True

# === Lancer Barrnap sur tous les bins ===
for binner, info in binners_info.items():
    binner_dir = os.path.join(args.binning_root, binner)
    if not os.path.exists(binner_dir):
        continue

    for sample in os.listdir(binner_dir):
        sample_dir = os.path.join(binner_dir, sample, info['pattern']) if info['pattern'] else os.path.join(binner_dir, sample)
        if not os.path.isdir(sample_dir):
            continue

        # Créer le dossier de sortie pour ce binner-sample
        output_dir = os.path.join(args.output_root, binner, sample)
        os.makedirs(output_dir, exist_ok=True)

        for bin_file in os.listdir(sample_dir):
            if not bin_file.endswith(info['ext']):
                continue

            bin_path = os.path.join(sample_dir, bin_file)
            fasta_tmp = os.path.join(output_dir, bin_file.replace(info['ext'], '.fasta'))
            gff_output = fasta_tmp.replace('.fasta', '.gff')

            # Extraction ou copie du fasta
            if bin_path.endswith('.gz'):
                with gzip.open(bin_path, 'rb') as f_in, open(fasta_tmp, 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)
            else:
                shutil.copy(bin_path, fasta_tmp)

            # Barrnap
            run_barrnap_on_file(fasta_tmp, gff_output)

# === Correction du summary final à partir des GFF ===
print("\n🔍 Calcul du vrai résumé 16S à partir des fichiers GFF...")

binner_samples = []
for binner in os.listdir(args.output_root):
    binner_path = os.path.join(args.output_root, binner)
    if not os.path.isdir(binner_path):
        continue

    for sample in os.listdir(binner_path):
        sample_path = os.path.join(binner_path, sample)
        if not os.path.isdir(sample_path):
            continue

        gff_files = [f for f in os.listdir(sample_path) if f.endswith('.gff')]

        n_total = len(gff_files)
        n_with_16s = 0

        for gff_file in gff_files:
            gff_path = os.path.join(sample_path, gff_file)
            with open(gff_path, 'r') as f:
                for line in f:
                    if not line.startswith('#') and '16S' in line:
                        n_with_16s += 1
                        break  # un seul comptage par fichier

        percent = (n_with_16s / n_total * 100) if n_total else 0
        binner_samples.append({
            'binner': binner,
            'sample': sample,
            'n_bins_total': n_total,
            'n_bins_with_16s': n_with_16s,
            'percent_with_16s': round(percent, 2),
        })

# === Sauvegarde finale du summary ===
df_summary = pd.DataFrame(binner_samples)
summary_file = os.path.join(args.output_root, 'summary_16S_detection.tsv')
df_summary.to_csv(summary_file, sep='\t', index=False)

print(f"\n✅ Résumé 16S  sauvegardé dans : {summary_file}")
print(f"🗂️ Dossiers des résultats par binner/sample dans : {args.output_root}")
