#!/usr/bin/env python3

import os
import argparse
from Bio import SeqIO
import pandas as pd

PADDING = 100

def get_bin_number(file_name, binner, sample):
    if binner == "maxbin":
        return file_name.replace(sample + "_maxbin.", "").replace(".fasta", "")
    elif binner == "semibin":
        return file_name.replace("SemiBin_", "").replace(".fasta", "")
    elif binner == "vamb":
        return file_name.replace(".fasta", "")
    return file_name.replace(".fasta", "")

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-g', '--gff_root', required=True, help="Racine contenant les GFF et FASTA par {binner}/{sample}")
    parser.add_argument('-o', '--output_root', required=True, help="Dossier de sortie des séquences 16S extraites")
    args = parser.parse_args()

    os.makedirs(args.output_root, exist_ok=True)
    all_hits = []
    failures = []

    for binner in os.listdir(args.gff_root):
        binner_path = os.path.join(args.gff_root, binner)
        if not os.path.isdir(binner_path):
            continue

        for sample in os.listdir(binner_path):
            sample_path = os.path.join(binner_path, sample)
            if not os.path.isdir(sample_path):
                continue

            output_dir = os.path.join(args.output_root, binner, sample)
            os.makedirs(output_dir, exist_ok=True)

            fasta_files = [f for f in os.listdir(sample_path) if f.endswith('.fasta')]
            for fasta_file in fasta_files:
                fasta_path = os.path.join(sample_path, fasta_file)
                bin_number = get_bin_number(fasta_file, binner, sample)
                gff_file = fasta_file.replace('.fasta', '.gff')
                gff_path = os.path.join(sample_path, gff_file)

                if not os.path.exists(gff_path):
                    failures.append({'binner': binner, 'sample': sample, 'bin': bin_number, 'reason': 'GFF not found'})
                    continue

                try:
                    records = {r.id: r for r in SeqIO.parse(fasta_path, 'fasta')}
                except Exception as e:
                    failures.append({'binner': binner, 'sample': sample, 'bin': bin_number, 'reason': f'FASTA parse error: {e}'})
                    continue

                seen_ids = set()
                output_fasta_path = os.path.join(output_dir, f"{binner}_{sample}_{bin_number}_16S_extracted.fasta")
                found = False
                lines_to_write = []

                with open(gff_path) as gff:
                    for line in gff:
                        if line.startswith('#') or not line.strip():
                            continue
                        fields = line.strip().split('\t')
                        if len(fields) < 9 or '16S' not in fields[8]:
                            continue

                        contig, start, end, strand = fields[0], int(fields[3]), int(fields[4]), fields[6]
                        if contig not in records:
                            failures.append({'binner': binner, 'sample': sample, 'bin': bin_number, 'reason': f'Contig {contig} not found'})
                            continue

                        uid = f"{binner}_{sample}_{bin_number}|{contig}_{start}_{end}_{strand}"
                        if uid in seen_ids:
                            continue
                        seen_ids.add(uid)

                        seq = records[contig].seq[max(0, start - 1 - PADDING): end + PADDING]
                        if strand == '-':
                            seq = seq.reverse_complement()

                        lines_to_write.append(f">{uid}\n{seq}\n")
                        all_hits.append({
                            'binner': binner,
                            'sample': sample,
                            'bin': bin_number,
                            'contig': contig,
                            'start': start,
                            'end': end,
                            'strand': strand
                        })
                        found = True

                if found:
                    with open(output_fasta_path, 'w') as fout:
                        fout.writelines(lines_to_write)
                else:
                    failures.append({'binner': binner, 'sample': sample, 'bin': bin_number, 'reason': 'no 16S found'})

    # === Sauvegarde des résultats ===
    if all_hits:
        pd.DataFrame(all_hits).to_csv(os.path.join(args.output_root, "16S_hits.tsv"), sep='\t', index=False)

    if failures:
        fail_dir = os.path.join(args.output_root, "failures")
        os.makedirs(fail_dir, exist_ok=True)
        pd.DataFrame(failures).to_csv(os.path.join(fail_dir, "failures.tsv"), sep='\t', index=False)
        print(f"⚠️ Fichiers échoués sauvegardés dans : {fail_dir}/failures.tsv")

    print("✅ Extraction 16S terminée.")

if __name__ == "__main__":
    main()
