#!/usr/bin/env python

import os
import subprocess
import argparse
import gzip
import shutil
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument('--binning_root', required=True)
parser.add_argument('--output_root', required=True)
args = parser.parse_args()

binners_info = {
    'semibin': {'pattern': 'output_bins', 'ext': '.fa.gz'},
    'maxbin': {'pattern': '', 'ext': '.fasta'},
    'vamb': {'pattern': 'bins', 'ext': '.fna'},
}

summary = []

def run_barrnap_on_file(fasta_in, gff_out):
    cmd = f"barrnap --kingdom bac {fasta_in} > {gff_out}"
    result = subprocess.run(cmd, shell=True)
    return result.returncode == 0

for binner, info in binners_info.items():
    binner_dir = os.path.join(args.binning_root, binner)
    if not os.path.exists(binner_dir):
        continue
    for sample in os.listdir(binner_dir):
        sample_dir = os.path.join(binner_dir, sample, info['pattern']) if info['pattern'] else os.path.join(binner_dir, sample)
        if not os.path.isdir(sample_dir):
            continue

        output_dir = os.path.join(args.output_root, binner, sample)
        os.makedirs(output_dir, exist_ok=True)

        n_bins_total, n_bins_with_16s = 0, 0
        for bin_file in os.listdir(sample_dir):
            if not bin_file.endswith(info['ext']):
                continue

            bin_path = os.path.join(sample_dir, bin_file)
            fasta_tmp = os.path.join(output_dir, bin_file.replace(info['ext'], '.fasta'))
            gff_output = fasta_tmp.replace('.fasta', '.gff')

            if bin_path.endswith('.gz'):
                with gzip.open(bin_path, 'rb') as f_in, open(fasta_tmp, 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)
            else:
                shutil.copy(bin_path, fasta_tmp)

            found = run_barrnap_on_file(fasta_tmp, gff_output)

            n_bins_total += 1
            if found:
                n_bins_with_16s += 1

        percent = (n_bins_with_16s / n_bins_total * 100) if n_bins_total else 0
        summary.append({
            'binner': binner, 'sample': sample,
            'n_bins_total': n_bins_total, 'n_bins_with_16s': n_bins_with_16s,
            'percent_with_16s': round(percent, 2),
        })

df_summary = pd.DataFrame(summary)
os.makedirs(args.output_root, exist_ok=True)
df_summary.to_csv(os.path.join(args.output_root, 'summary_16S_detection.tsv'), sep='\t', index=False)

print("\n✅ Analyse 16S terminée proprement pour tous les bins.")
