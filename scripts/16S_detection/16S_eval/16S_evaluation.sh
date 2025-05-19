#!/bin/bash
#SBATCH --job-name=16S_eval_array
#SBATCH --array=0-6
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=logs/16S_eval_1sample_1binner_focus_%A_%a.out
#SBATCH --error=logs/16S_eval_1sample_1binner_focus_%A_%a.err

echo "🚀 Début du job array $SLURM_ARRAY_TASK_ID"

# === Activation de l’environnement ===
source activate mahabio_env

# === Lire la ligne correspondant à ce job array ===
BINNER_SAMPLE=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" /shared/home/asandri/MAHABIO_analysis/scripts/16S_detection/old_scripts/samples_list.txt)

BINNER=$(echo "$BINNER_SAMPLE" | cut -d'/' -f1)
SAMPLE=$(echo "$BINNER_SAMPLE" | cut -d'/' -f2)

FASTA_DIR="/shared/home/asandri/MAHABIO_analysis/results/binning/Affiliation/16S_extracted/${BINNER}/${SAMPLE}"
EVAL_DIR="/shared/home/asandri/MAHABIO_analysis/results/binning/Affiliation/16S_evaluation"
ALIGN_DIR="$EVAL_DIR/aligned/${SAMPLE}/${BINNER}"
PLOT_DIR="$EVAL_DIR/plots/${SAMPLE}/${BINNER}"
IDENTITY_DIR="$EVAL_DIR/identity/${SAMPLE}/${BINNER}"
DIVERG_DIR="$EVAL_DIR/divergence_logs/${SAMPLE}/${BINNER}/"
mkdir -p "$ALIGN_DIR" "$PLOT_DIR" "$IDENTITY_DIR" "$DIVERG_DIR"

# === Étape 1 : filtrer les bins divergents ===
echo "🔍 [${BINNER}/${SAMPLE}] Filtrage des bins divergents"
python /shared/home/asandri/MAHABIO_analysis/scripts/16S_detection/16S_eval/filter_divergent_bins.py \
  --fasta_dir "$FASTA_DIR" \
  --threshold 0.97 > "$DIVERG_DIR/${BINNER}_${SAMPLE}.log"

# === Étape 2 : aligner les séquences par bin ===
echo "🧬 [${BINNER}/${SAMPLE}] Alignement par bin"
python /shared/home/asandri/MAHABIO_analysis/scripts/16S_detection/16S_eval/align_16S_per_bin.py \
  --fasta_dir "$FASTA_DIR" \
  --output_dir "$ALIGN_DIR"

# === Étape 3 : visualisation des alignements ===
echo "🖼️ [${BINNER}/${SAMPLE}] Génération des visualisations"
python /shared/home/asandri/MAHABIO_analysis/scripts/16S_detection/16S_eval/plot_16S_alignment.py \
  -i "$ALIGN_DIR" -o "$PLOT_DIR"

# === Étape 4 : matrice d’identité ===
echo "📊 [${BINNER}/${SAMPLE}] Calcul des matrices d'identité"
python /shared/home/asandri/MAHABIO_analysis/scripts/16S_detection/16S_eval/compute_identity_matrix.py \
  -i "$ALIGN_DIR" -o "$IDENTITY_DIR"

echo "✅ [${BINNER}/${SAMPLE}] Terminé"
