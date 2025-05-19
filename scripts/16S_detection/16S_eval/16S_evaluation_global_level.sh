#!/bin/bash
#SBATCH --job-name=16S_eval_all
#SBATCH --output=logs/16S_eval_all_%A_%a.out
#SBATCH --error=logs/16S_eval_all_%A_%a.err
#SBATCH --time=03:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --array=0-2

# === Paramètres ===
samples=("C_contigs_more_than_300bp" "Fk_contigs_more_than_300bp" "Sj_contigs_more_than_300bp")
SAMPLE=${samples[$SLURM_ARRAY_TASK_ID]}

# === Activation de l’environnement ===
source activate mahabio_env

# === Chemins ===
BASE_DIR="/shared/home/asandri/MAHABIO_analysis/results/binning/Affiliation"
FASTA_DIR="$BASE_DIR/16S_extracted"
EVAL_DIR="$BASE_DIR/16S_evaluation"
ALIGN_DIR="$EVAL_DIR/aligned/${SAMPLE}"
PLOT_DIR="$EVAL_DIR/plots/${SAMPLE}/"
IDENTITY_DIR="$EVAL_DIR/identity/"
DIVERG_DIR="$EVAL_DIR/divergence_logs/"
mkdir -p "$ALIGN_DIR" "$PLOT_DIR" "$IDENTITY_DIR" "$DIVERG_DIR"

# === Étape 1 : concaténation des 16S tous binners pour le sample ===
python align_16S_per_sample_all_bin_then_all_binners.py \
  --fasta_dir "$FASTA_DIR" \
  --sample "$SAMPLE" \
  --output_dir "$ALIGN_DIR"

# === Étape 2 : filtrage divergence ===
python filter_divergent_bins.py \
  --fasta_dir "$ALIGN_DIR" \
  --threshold 0.97 > "$DIVERG_DIR/${SAMPLE}_ALL_BINS.log"

# === Étape 3 : visualisation alignements ===
python plot_16S_alignment.py \
  -i "$ALIGN_DIR" -o "$PLOT_DIR"

# === Étape 4 : matrices d’identité ===
python compute_identity_matrix.py \
  -i "$ALIGN_DIR" -o "$IDENTITY_DIR"
#visualise heatmap for matrices 
python heatmap_identity.py


# === Étape 5 : affiliation taxonomique (préparation BLASTN) ===
echo "> Préparation de l'affiliation taxonomique... (blastn ou GraftM à intégrer ensuite)"

# À compléter : blastn contre SILVA ou GraftM sur les séquences concaténées
# blastn -query ${ALIGN_DIR}/${SAMPLE}_ALL_BINS.fasta -db SILVA -out ...
# graftM graft --input ...

echo "✅ Pipeline 16S global terminé pour $SAMPLE"
