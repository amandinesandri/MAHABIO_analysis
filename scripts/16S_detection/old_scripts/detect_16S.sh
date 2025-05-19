#!/bin/bash
#SBATCH --job-name=detect_16S_barrnap
#SBATCH --output=logs/detect_16S_barrnap_%j.out
#SBATCH --error=logs/detect_16S_barrnap_%j.err
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G

# === Activation de l'environnement ===
echo "🚀 Activation de l'environnement conda detect16S_env"
source activate mahabio_env

# === Création des dossiers ===
echo "📂 Préparation des répertoires de sortie"
mkdir -p /shared/home/asandri/MAHABIO_analysis/results/binning/Affiliation/16S_detection/

# === Lancement du script principal ===
echo "🔍 Lancement de la détection 16S avec Barrnap"
python detect16S.py \
  --binning_root /shared/home/asandri/MAHABIO_analysis/results/binning \
  --output_root /shared/home/asandri/MAHABIO_analysis/results/binning/Affiliation/16S_detection

# === Fin ===
echo "✅ Détection 16S terminée et outputs correctement organisés !"

# python /shared/home/asandri/MAHABIO_analysis/scripts/Binning/recalcul_summary_16s_from_gff.py \
#   --gff_root /shared/home/asandri/MAHABIO_analysis/results/binning/Affiliation/16S_detection \
#   --output_summary /shared/home/asandri/MAHABIO_analysis/results/binning/Affiliation/16S_detection/summary_16S_detection_corrected.tsv
