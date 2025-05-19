#!/bin/bash
#SBATCH --job-name=hmmsearch_cmuA
#SBATCH --output=logs/hmmsearch/hmmsearch_cmuA_%A_%a.out
#SBATCH --error=logs/hmmsearch/hmmsearch_cmuA_%A_%a.err
#SBATCH --time=01:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=2
#SBATCH --array=0-389

source activate mahabio_env

# === Variables
BASE="/shared/home/asandri/MAHABIO_analysis"
BINS_DIR="$BASE/results/binning/checkm"
HMM="$BASE/results/hmm/cmuA.hmm"
OUT_BASE="$BASE/results/binning/cmuA_hmmsearch"
mkdir -p "$OUT_BASE"

# === Récupérer tous les genes.faa des bins
FILES=($(find "$BINS_DIR" -type f -name "genes.faa"))

# === Sélectionner le fichier correspondant
faa="${FILES[$SLURM_ARRAY_TASK_ID]}"

# === Extraire les infos
rel_path=${faa#$BINS_DIR/}
binner=$(echo "$rel_path" | cut -d'/' -f1)
sample=$(echo "$rel_path" | cut -d'/' -f2)
bin_folder=$(basename "$(dirname "$faa")")  # ATTENTION: c'est le dossier contenant genes.faa

# === Définir bin_number (toujours le numéro pur)
if [[ "$binner" == "maxbin" ]]; then
    bin_number=$(echo "$bin_folder" | sed -E 's/.*\.([0-9]+)$/\1/')
elif [[ "$binner" == "semibin" ]]; then
    bin_number=$(echo "$bin_folder" | sed -E 's/.*SemiBin_([0-9]+)$/\1/')
elif [[ "$binner" == "vamb" ]]; then
    bin_number="$bin_folder"
else
    bin_number="$bin_folder"
fi

# === Construire ref_bin de manière cohérente
ref_bin="${binner}_${sample}_${bin_number}"

# === Créer le dossier de sortie
OUT_DIR="$OUT_BASE/$binner/$sample"
mkdir -p "$OUT_DIR"

# === Lancer hmmsearch
echo "🔎 HMMsearch sur bin: $ref_bin"
hmmsearch --tblout "$OUT_DIR/${ref_bin}_cmuA.tbl" "$HMM" "$faa"

echo "✅ Fini pour bin: $ref_bin"
