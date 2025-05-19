#!/bin/bash
#SBATCH --job-name=master_itol_pipeline
#SBATCH --output=logs/master_itol_pipeline_%j.out
#SBATCH --error=logs/master_itol_pipeline_%j.err
#SBATCH --time=00:10:00
#SBATCH --mem=2G
#SBATCH --cpus-per-task=1




source ~/miniforge3/bin/activate mahabio_env
echo "🚀 Lancement du pipeline iTOL complet..."

# Reset nettoyage avant relancer master_pipeline_itol.sh

# Activer l'environnement

# Définir la base
BASE="/shared/home/asandri/MAHABIO_analysis"
ALIGN_DIR="$BASE/results/binning/cmuA_hmmsearch/alignment_phylo"
HITS_DIR="$BASE/results/binning/cmuA_hmmsearch/cmuA_hits_sequences"


# Supprimer les anciens fichiers générés
rm -f "$ALIGN_DIR/all_cmuA_hits_raw.faa" \
      "$ALIGN_DIR/all_cmuA_hits_unique.faa" \
      "$ALIGN_DIR/all_cmuA_hits_aligned.faa" \
      "$ALIGN_DIR/all_cmuA_hits_tree.nwk" \
      "$ALIGN_DIR/all_cmuA_hits_tree_rerooted.nwk" \
      "$ALIGN_DIR/cmuA_mapping.tsv" \
      "$ALIGN_DIR/itol_cmuA_metadata_final_fixed.txt" \
      "$ALIGN_DIR/itol_functional_groups.txt" \
      "$ALIGN_DIR/itol_functional_groups_final_fixed.txt"

# Vider le dossier des séquences extraites
rm -f "$HITS_DIR"/*.faa

echo "✅ Tous les anciens fichiers générés ont été supprimés proprement."
echo "🚀 Prêt à relancer master_pipeline_itol.sh !"

# Lancement séquentiel
sbatch hmmsearch_cmua.sh
sleep 2

sbatch parse_cmua_hits.sh
sleep 2

sbatch extract_hits.sh
sleep 2

sbatch Alignement.sh
sleep 2

sbatch reroot_tree.sh
sleep 2

sbatch generate_cmuA_mapping.sh
sleep 2

sbatch generate_itol_mapping_final.sh
sleep 2

sbatch generate_itol_functional_groups_final.sh
sleep 2

echo "🎯 Tous les jobs sont soumis ! Suivi via squeue -u \$USER"

 Dossiers
BASE="/shared/home/asandroutput_pathi/MAHABIO_analysis"
ALIGN_DIR="$BASE/results/binning/cmuA_hmmsearch/alignment_phylo"
HITS_DIR="$BASE/results/binning/cmuA_hmmsearch/cmuA_hits_sequences"

# Fichiers attendus
expected_files=(
    "$ALIGN_DIR/all_cmuA_hits_raw.faa"
    "$ALIGN_DIR/all_cmuA_hits_unique.faa"
    "$ALIGN_DIR/all_cmuA_hits_aligned.faa"
    "$ALIGN_DIR/all_cmuA_hits_tree.nwk"
    "$ALIGN_DIR/all_cmuA_hits_tree_rerooted.nwk"
    "$ALIGN_DIR/cmuA_mapping.tsv"
    "$ALIGN_DIR/itol_cmuA_metadata_final_fixed.txt"
    "$ALIGN_DIR/itol_functional_groups.txt"
    "$ALIGN_DIR/itol_functional_groups_final_fixed.txt"
)

# Check
echo "🔎 Vérification des fichiers générés dans : $ALIGN_DIR"
for file in "${expected_files[@]}"; do
    if [[ -s "$file" ]]; then
        echo "✅ $file existe et n'est pas vide."
    elif [[ -e "$file" ]]; then
        echo "⚠️ $file existe mais est vide !"
    else
        echo "❌ $file est manquant !"
    fi
done

# Vérifier aussi s'il y a bien des séquences extraites
echo ""
echo "🔎 Vérification des séquences extraites dans : $HITS_DIR"
n_faa=$(ls "$HITS_DIR"/*.faa 2>/dev/null | wc -l)

if [[ "$n_faa" -gt 0 ]]; then
    echo "✅ $n_faa fichiers hits extraits présents."
else
    echo "❌ Aucun fichier de hits extrait trouvé !"
fi

echo "🎯 Vérification terminée."