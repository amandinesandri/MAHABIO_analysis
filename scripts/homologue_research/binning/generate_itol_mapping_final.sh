# generate_itol_mapping_final.sh (version corrigée et enrichie)

#!/bin/bash
#SBATCH --job-name=itol_mapping_final
#SBATCH --output=logs/itol_mapping_final_%j.out
#SBATCH --error=logs/itol_mapping_final_%j.err
#SBATCH --time=00:20:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1

source ~/miniforge3/bin/activate mahabio_env

# Variables
BASE="/shared/home/asandri/MAHABIO_analysis"
ALIGN_FASTA="$BASE/results/binning/cmuA_hmmsearch/alignment_phylo/all_cmuA_hits_aligned.faa"
MAPPING_FILE="$BASE/results/binning/cmuA_hmmsearch/alignment_phylo/cmuA_mapping_enriched.tsv"
ITOL_MAPPING_OUT="$BASE/results/binning/cmuA_hmmsearch/alignment_phylo/itol_cmuA_metadata.txt"

# Script Python inline
python - <<EOF
import pandas as pd

# Charger les données enrichies
mapping_df = pd.read_csv("$MAPPING_FILE", sep="\t")

# Lire les IDs présents dans l'alignement
with open("$ALIGN_FASTA") as f:
    seq_ids = [line[1:].strip() for line in f if line.startswith(">")]

# Couleurs associées aux binners
color_map = {
    "maxbin": "#1f77b4",
    "semibin": "#2ca02c",
    "vamb": "#d62728"
}

# Créer le fichier iTOL enrichi
with open("$ITOL_MAPPING_OUT", "w") as out:
    out.write("DATASET_COLORSTRIP\n")
    out.write("SEPARATOR TAB\n")
    out.write("DATASET_LABEL\tBinner_and_Properties\n")
    out.write("COLOR\t#ff0000\n")
    out.write("DATA\n")

    for seq_id in seq_ids:
        row = mapping_df[mapping_df['seq_id'] == seq_id]
        if not row.empty:
            binner = row.iloc[0]['binner']
            strand = row.iloc[0].get('strand', 'NA')
            partial = row.iloc[0].get('partial', 'NA')
            gc_content = row.iloc[0].get('gc_content', 'NA')
            color = color_map.get(binner.lower(), "#7f7f7f")

            label = f"{binner}_strand:{strand}_partial:{partial}_GC:{gc_content}"
            out.write(f"{seq_id}\t{color}\t{label}\n")
        else:
            print(f"Warning: {seq_id} not found in mapping table!")
EOF

# Message final
echo "\ud83d\udc9a iTOL mapping enrichi généré : \$ITOL_MAPPING_OUT"
