# generate_itol_functional_groups_final.sh

#!/bin/bash
#SBATCH --job-name=itol_functional_final
#SBATCH --output=logs/itol_functional_final_%j.out
#SBATCH --error=logs/itol_functional_final_%j.err
#SBATCH --time=00:20:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1

source ~/miniforge3/bin/activate mahabio_env

# Variables
BASE="/shared/home/asandri/MAHABIO_analysis"
ALIGN_FASTA="$BASE/results/binning/cmuA_hmmsearch/alignment_phylo/all_cmuA_hits_aligned.faa"
FUNCTIONAL_RAW="$BASE/results/binning/cmuA_hmmsearch/alignment_phylo/itol_functional_groups.txt"
ITOL_FUNCTIONAL_OUT="$BASE/results/binning/cmuA_hmmsearch/alignment_phylo/itol_functional_groups_final_fixed.txt"

# Script Python inline
python - <<EOF
import pandas as pd

# Lire séquences présentes
seq_ids = []
with open("$ALIGN_FASTA") as f:
    for line in f:
        if line.startswith(">"):
            seq_ids.append(line[1:].strip())

# Charger le fichier fonctionnel
func_df = pd.read_csv("$FUNCTIONAL_RAW", sep="\t", comment="#", names=["seq_id", "color"], skiprows=5)
seq_to_color = dict(zip(func_df['seq_id'], func_df['color']))

# Générer iTOL
with open("$ITOL_FUNCTIONAL_OUT", "w") as out:
    out.write("DATASET_COLORSTRIP\n")
    out.write("SEPARATOR TAB\n")
    out.write("DATASET_LABEL\tFunctional Groups\n")
    out.write("COLOR\t#0000ff\n")
    out.write("DATA\n")
    for seq_id in seq_ids:
        color = seq_to_color.get(seq_id)
        if color:
            out.write(f"{seq_id}\t{color}\n")
        else:
            print(f"Warning: {seq_id} not found in functional groups!")
EOF

echo "✅ iTOL Functional Groups fichier généré : $ITOL_FUNCTIONAL_OUT"
