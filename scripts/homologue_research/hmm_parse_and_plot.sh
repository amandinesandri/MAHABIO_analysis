#!/bin/bash
#SBATCH --job-name=parse_hmm
#SBATCH --output=logs/parse_hmm_%j.out
#SBATCH --error=logs/parse_hmm_%j.err
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G

# === 1. Activation de l’environnement ===
source activate mahabio_env

# === 2. Définition des chemins ===
DOMTBL_DIR="/shared/home/asandri/MAHABIO_analysis/results/hmm"
SUMMARY_DIR="/shared/home/asandri/MAHABIO_analysis/results/hmm/summary"
COMBINED="${SUMMARY_DIR}/combined_hits_summary.tsv"
PLOT_OUT="/shared/home/asandri/MAHABIO_analysis/results/hmm/plots"
CHECK_SCRIPT="/shared/home/asandri/MAHABIO_analysis/scripts/homologue_research/check_tsv_integrity.sh"

mkdir -p "$SUMMARY_DIR" "$PLOT_OUT" logs

# === 3. Parsing + nettoyage de chaque .domtbl ===
for domtbl in "$DOMTBL_DIR"/*_hits.domtbl; do
    base=$(basename "$domtbl" _hits.domtbl)
    out_tsv="${SUMMARY_DIR}/${base}_summary.tsv"
    echo "📄 Parsing $domtbl → $out_tsv"
    
    python /shared/home/asandri/MAHABIO_analysis/scripts/homologue_research/hmm_parse.py \
        "$domtbl" "$out_tsv"
    
    # === Nettoyage : on garde uniquement les lignes avec >1 champ ===
    awk -F'\t' 'NF>1' "$out_tsv" > "${out_tsv}.clean"
    mv "${out_tsv}.clean" "$out_tsv"
done

# === 4. Concaténation robuste en un seul fichier ===
echo -e "#sample\t$(head -1 ${SUMMARY_DIR}/*_summary.tsv | head -n1)" > "$COMBINED"
for f in ${SUMMARY_DIR}/*_summary.tsv; do
    sample=$(basename "$f" _summary.tsv)
    
    # Skip les fichiers sans données
    if [[ $(wc -l < "$f") -lt 2 ]]; then
        echo "⚠️ Skipping $sample: no valid data found."
        continue
    fi
    
    tail -n +2 "$f" | awk -v sample_name="$sample" 'BEGIN{FS=OFS="\t"}{print sample_name, $0}' >> "$COMBINED"

done
echo "✅ Fichier combiné : $COMBINED"

# === 5. Vérification d'intégrité TSV ===
echo "🔎 Vérification du fichier TSV combiné..."
bash "$CHECK_SCRIPT" "$COMBINED"
STATUS=$?

if [[ "$STATUS" -ne 0 ]]; then
    echo "❌ Erreur : fichier TSV non conforme, arrêt du pipeline."
    exit 1
fi

# === 6. Création du plot matplotlib ===
echo "📊 Génération du plot des hits par échantillon…"
python - << EOF
import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv("${COMBINED}", sep="\t")
counts = df['sample'].value_counts().sort_index()

plt.figure(figsize=(6,4))
counts.plot(kind='bar')
plt.title("Nombre de hits HMM par échantillon")
plt.xlabel("Échantillon")
plt.ylabel("Nombre de hits")
plt.tight_layout()
plt.savefig("${PLOT_OUT}/hmm_hits_per_sample.png")
EOF

# === 7. Final ===
echo "🎉 Tout est terminé ! Résultats disponibles :"
echo "  - Résumés individuels : $SUMMARY_DIR/"
echo "  - Fichier combiné : $COMBINED"
echo "  - Figure : $PLOT_OUT/hmm_hits_per_sample.png"
