#!/bin/bash
#SBATCH --job-name=extract_cmuA_hits
#SBATCH --output=logs/extract_cmuA_hits_%j.out
#SBATCH --error=logs/extract_cmuA_hits_%j.err
#SBATCH --time=01:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=2

# Activer l'environnement
source activate mahabio_env

# Variables
BASE="/shared/home/asandri/MAHABIO_analysis"
HITS_TSV="$BASE/results/binning/cmuA_hmmsearch/cmuA_hits_summary_1e-5_50.tsv"
OUT_DIR="$BASE/results/binning/cmuA_hmmsearch/cmuA_hits_sequences"
BINS_DIR="$BASE/results/binning/checkm"

mkdir -p "$OUT_DIR"

# Lire et traiter chaque ligne du fichier des hits
tail -n +2 "$HITS_TSV" | while IFS=$'\t' read -r binner sample bin_number ref_bin target_name query_name evalue score; do

    # Définir faa_path en fonction du binner
    if [[ "$binner" == "maxbin" ]]; then
        faa_path="$BINS_DIR/maxbin/$sample/bins/${sample}_maxbin.${bin_number}/genes.faa"
    elif [[ "$binner" == "semibin" ]]; then
        faa_path="$BINS_DIR/semibin/$sample/bins/SemiBin_${bin_number}/genes.faa"
    elif [[ "$binner" == "vamb" ]]; then
        faa_path="$BINS_DIR/vamb/$sample/bins/${bin_number}/genes.faa"
    else
        faa_path=""
    fi

    output_path="$OUT_DIR/${ref_bin}_cmuA_hits.faa"

    if [[ -f "$faa_path" ]]; then
        seqkit grep -p "$target_name" "$faa_path" >> "$output_path"
        echo "✅ Séquence $target_name extraite dans $output_path"
    else
        echo "⚠️ Fichier $faa_path introuvable !"
    fi

done

echo "🎯 Extraction terminée. Résultats dans : $OUT_DIR"
