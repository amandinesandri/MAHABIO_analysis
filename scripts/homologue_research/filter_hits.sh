#!/bin/bash
#SBATCH --job-name=filter_hits_per_sample
#SBATCH --output=logs/filter_hits_per_sample_%j.out
#SBATCH --error=logs/filter_hits_per_sample_%j.err
#SBATCH --time=00:20:00
#SBATCH --mem=4G

source activate mahabio_env

# Variables
BASE="/shared/home/asandri/MAHABIO_analysis"
SUMMARY_DIR="$BASE/results/hmm/summary"
FILTERED_DIR="$BASE/results/hmm/filtered_hits"
mkdir -p "$FILTERED_DIR"

# Boucle sur chaque fichier summary
for summary in "$SUMMARY_DIR"/*_summary.tsv; do
    sample=$(basename "$summary" _summary.tsv)
    output="$FILTERED_DIR/${sample}_filtered_hits.tsv"

    echo "Filtrage de $sample"

    python - << EOF
import pandas as pd
df = pd.read_csv("${summary}", sep="\t")
filtered = df[
    (df['evalue'] < 1e-5) &
    (df['score'] > 50) &
    (df['coverage'] > 0.30)
]
filtered.to_csv("${output}", sep="\t", index=False)
print(f"✅ {len(filtered)} hits retenus pour ${sample}")
EOF
done
