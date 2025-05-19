#!/bin/bash
#SBATCH --job-name=generate_cmuA_mapping_enriched
#SBATCH --output=logs/generate_cmuA_mapping_enriched_%j.out
#SBATCH --error=logs/generate_cmuA_mapping_enriched_%j.err
#SBATCH --time=00:15:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1

source ~/miniforge3/bin/activate mahabio_env

# Variables
BASE="/shared/home/asandri/MAHABIO_analysis"
HITS_TSV="$BASE/results/binning/cmuA_hmmsearch/cmuA_hits_summary_1e-5_50.tsv"
HITS_DIR="$BASE/results/binning/cmuA_hmmsearch/cmuA_hits_sequences"
MAPPING_OUT="$BASE/results/binning/cmuA_hmmsearch/alignment_phylo/cmuA_mapping_enriched.tsv"

#python generate_cmuA_mapping_enriched.py "$HITS_TSV" "$HITS_DIR" "$MAPPING_OUT"


# Script Python inline
python - <<EOF
import os
import pandas as pd

# Charger hits cmuA
hits = pd.read_csv("$HITS_TSV", sep="\t")
hits['seq_id'] = ['seq_' + str(i+1) for i in range(len(hits))]

# Fonction pour extraire les infos du header fasta
def extract_header_info(header):
    parts = header.split("#")
    if len(parts) >= 5:
        try:
            start = int(parts[1].strip())
            end = int(parts[2].strip())
            strand = parts[3].strip()
            attributes = parts[4]
            partial = None
            gc_content = None
            for item in attributes.split(";"):
                if item.startswith("partial="):
                    partial = item.split("=")[1]
                if item.startswith("gc_cont="):
                    gc_content = float(item.split("=")[1])
            return start, end, strand, partial, gc_content
        except Exception:
            return None, None, None, None, None
    else:
        return None, None, None, None, None

# Parcourir et enrichir
start_list = []
end_list = []
strand_list = []
partial_list = []
gc_content_list = []

for idx, row in hits.iterrows():
    ref_bin = row['ref_bin']
    target_name = row['target_name']
    faa_path = os.path.join("$HITS_DIR", f"{ref_bin}_cmuA_hits.faa")

    if os.path.isfile(faa_path):
        found = False
        with open(faa_path) as f:
            for line in f:
                if line.startswith(">"):
                    header_core = line[1:].split(" ")[0]
                    if header_core == target_name:
                        start, end, strand, partial, gc_content = extract_header_info(line[1:].strip())
                        start_list.append(start)
                        end_list.append(end)
                        strand_list.append(strand)
                        partial_list.append(partial)
                        gc_content_list.append(gc_content)
                        found = True
                        break
        if not found:
            start_list.append(None)
            end_list.append(None)
            strand_list.append(None)
            partial_list.append(None)
            gc_content_list.append(None)
    else:
        start_list.append(None)
        end_list.append(None)
        strand_list.append(None)
        partial_list.append(None)
        gc_content_list.append(None)

# Ajouter les colonnes
hits['start'] = start_list
hits['end'] = end_list
hits['strand'] = strand_list
hits['partial'] = partial_list
hits['gc_content'] = gc_content_list

# Sauver
hits.to_csv("$MAPPING_OUT", sep="\t", index=False)
print(f"\u2705 Mapping enrichi généré : {MAPPING_OUT}")
EOF

echo "\n\ud83c\udf89 Mapping enrichi OK dans : \$MAPPING_OUT"
