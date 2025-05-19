#!/bin/bash
#SBATCH --job-name=CH3Cl_array
#SBATCH --output=logs/CH3Cl_array_%A_%a.out
#SBATCH --error=logs/CH3Cl_array_%A_%a.err
#SBATCH --time=01:00:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1

source activate mahabio_env
# Fichier TSV passé en argument
BIN_FILE="$1"

# Récupération de la ligne correspondant à l’index SLURM_ARRAY_TASK_ID (ligne + entête)
LINE=$(awk -v line_num=$((SLURM_ARRAY_TASK_ID + 2)) 'NR == line_num' "$BIN_FILE")

# Extraction des champs tabulaires
IFS=$'\t' read -r ref_bin binner sample bin_number filename size_bp size_mbp size_class <<< "$LINE"

# Chemin vers le dossier du bin (on utilise filename !)
BIN_PATH="/shared/home/asandri/MAHABIO_analysis/results/binning/checkm/${binner}/${sample}/bins/${filename}"
mkdir -p /shared/home/asandri/MAHABIO_analysis/results/binning/Annotation/
# Choix du fichier d'annotation
if [ -f "${BIN_PATH}/genes.gff" ]; then
    ANNOT_FILE="${BIN_PATH}/genes.gff"
elif [ -f "${BIN_PATH}/hmmer.analyze.txt" ]; then
    ANNOT_FILE="${BIN_PATH}/hmmer.analyze.txt"
else
    echo -e "${binner}\t${sample}\t${filename}\tAUCUN_FICHIER" >> /shared/home/asandri/MAHABIO_analysis/results/binning/Annotation/log_bins_sans_annotation.tsv
    exit 0
fi

# Script Python inline pour extraire les gènes d’intérêt
python3 - <<EOF
import re, os
import pandas as pd

target_patterns = [
    "cmuA", "cmuB", "hmtA", "mtaA", "mtaB", "mtaC",
    "chloromethane", "corrinoid", "tetrahydrofolate",
    "methyltransferase", "methyl-H4F", "methyl-H4MPT"
]
patterns = [re.compile(rf"\\b{p}\\b", re.IGNORECASE) for p in target_patterns]

bin_path = "${BIN_PATH}"
annot_file = "${ANNOT_FILE}"

hits = []
with open(annot_file, "r") as f:
    for line in f:
        for pattern in patterns:
            if pattern.search(line):
                hits.append({
                    "binner": "${binner}",
                    "sample": "${sample}",
                    "bin_id": "${filename}",
                    "source_file": os.path.basename(annot_file),
                    "line": line.strip()
                })
                break

if hits:
    df = pd.DataFrame(hits)
    outdir = f"/shared/home/asandri/MAHABIO_analysis/results/binning/Annotation/${binner}/${sample}/bins/${filename}"
    os.makedirs(outdir, exist_ok=True)
    outfile = os.path.join(outdir, "CH3Cl_candidate_genes.tsv")
    df.to_csv(outfile, sep="\\t", index=False)
    print(f"✔️ Résultats écrits dans {outfile}")
else:
    with open("/shared/home/asandri/MAHABIO_analysis/results/binning/Annotation/log_bins_sans_resultat.tsv", "a") as logf:
        logf.write(f"${binner}\\t${sample}\\t${filename}\\tAUCUN_RESULTAT\\n")
    print(f"❗ Aucun gène CH3Cl détecté dans ${filename}")
EOF
