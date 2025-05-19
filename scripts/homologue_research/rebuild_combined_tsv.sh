#!/bin/bash
#SBATCH --job-name=rebuild_combined
#SBATCH --output=logs/rebuild_combined_%j.out
#SBATCH --error=logs/rebuild_combined_%j.err
#SBATCH --time=00:10:00
#SBATCH --mem=2G

# === Variables ===
BASE="/shared/home/asandri/MAHABIO_analysis"
OUT_DIR="$BASE/results/hmm/summary"
COMBINED="$OUT_DIR/combined_hits_summary.tsv"

# === 1. Nettoyage
rm -f "$COMBINED"

# === 2. Prendre l'en-tête d'un seul fichier exemple (le premier trouvé)
first_file=$(ls "$OUT_DIR"/*_summary.tsv | head -n 1)

# === 3. Construire l'en-tête combiné : sample + colonnes originales
echo -e "sample\t$(head -n 1 "$first_file")" > "$COMBINED"

# === 4. Ajouter les données de tous les fichiers
for f in "$OUT_DIR"/*_summary.tsv; do
    sample_name=$(basename "$f" _summary.tsv)
    
    if [[ $(wc -l < "$f") -lt 2 ]]; then
        echo "⚠️ Skip $sample_name (fichier vide)"
        continue
    fi
    
    # Ajouter la colonne sample au début de chaque ligne
    tail -n +2 "$f" | awk -v sample="$sample_name" 'BEGIN{OFS="\t"}{print sample, $0}' >> "$COMBINED"
done

echo "✅ Fichier combiné créé proprement : $COMBINED"
