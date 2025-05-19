#!/usr/bin/env python3

import os
import argparse
from Bio import SeqIO
import subprocess
from collections import defaultdict


def concat_16S_sequences(fasta_root, sample):
    """
    Concatène toutes les séquences 16S extraites pour un sample donné, 
    tous bins et tous binners confondus.
    Retourne un dictionnaire de listes par binner et une liste globale.
    """
    concat_by_binner = defaultdict(list)
    concat_all = []

    for binner in os.listdir(fasta_root):
        binner_path = os.path.join(fasta_root, binner)
        print(f"⚠️ check binner path : {binner_path}")

        sample_path = os.path.join(binner_path, sample)
        print(f"⚠️ check sample path : {sample_path}")
        if not os.path.isdir(sample_path):
            continue

        for fasta_file in os.listdir(sample_path):
            print(f"⚠️ check fasta file : {fasta_file}")
            if not fasta_file.endswith("_16S_extracted.fasta"):
                continue

            fasta_path = os.path.join(sample_path, fasta_file)
            print(f"⚠️ check fasta path : {fasta_path}")
            print(f"pourquoi remettre le sample dans le path...")
            records = list(SeqIO.parse(fasta_path, "fasta"))

            for record in records:
                # Ajoute info provenance dans l'ID
                #record.id = f"{binner}|{sample}|{record.id}"
                record.description = ""
                print(f"⚠️ check record : {record}")
                concat_by_binner[binner].append(record)
                concat_all.append(record)
                print(f"⚠️ check binner : {binner}")
                print(f"⚠️ check concat by binner path : {concat_by_binner}")
                print(f"⚠️ check concat_all : {concat_all}")


    return concat_by_binner, concat_all


def write_fasta(records, output_path):
    with open(output_path, "w") as f:
        SeqIO.write(records, f, "fasta")
        print(f"⚠️ check record to fasta  : {output_path}")


def align_with_mafft(input_fasta, output_fasta):
    if len(list(SeqIO.parse(input_fasta, "fasta"))) < 2:
        print(f"⚠️ Pas assez de séquences dans {input_fasta} pour alignement.")
        return
    cmd = f"mafft --auto {input_fasta} > {output_fasta}"
    subprocess.run(cmd, shell=True)
    print(f"✅ Alignement enregistré : {output_fasta}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--fasta_dir", required=True, help="Répertoire racine des fichiers 16S extraits par binner/sample")
    parser.add_argument("--sample", required=True, help="Nom du sample à traiter (ex: C, Fk, Sj)")
    parser.add_argument("--output_dir", required=True, help="Répertoire de sortie pour les alignements MAFFT")
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    concat_by_binner, concat_all = concat_16S_sequences(args.fasta_dir, args.sample)

    # Alignement par sample-binner (tous bins confondus)
    for binner, records in concat_by_binner.items():
        if len(records) < 2:
            print(f"⚠️ Pas assez de séquences pour {args.sample}-{binner}.")
            continue
        binner_fasta = os.path.join(args.output_dir, f"{args.sample}_{binner}_ALL_BINS_16S.fasta")
        aligned_fasta = binner_fasta.replace(".fasta", "_aligned.fasta")
        write_fasta(records, binner_fasta)
        align_with_mafft(binner_fasta, aligned_fasta)

    # Alignement global pour le sample (tous binners confondus)
    if len(concat_all) >= 2:
        all_fasta = os.path.join(args.output_dir, f"{args.sample}_ALL_BINS_ALL_BINNERS_16S.fasta")
        aligned_all = all_fasta.replace(".fasta", "_aligned.fasta")
        write_fasta(concat_all, all_fasta)
        align_with_mafft(all_fasta, aligned_all)
    else:
        print(f"⚠️ Moins de 2 séquences pour {args.sample} au total : pas d'alignement global.")


if __name__ == "__main__":
    main()
