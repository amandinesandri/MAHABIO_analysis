#!/bin/bash
#SBATCH --job-name=checkm_per_sample
#SBATCH --output=logs/checkm_per_sample_%j.out
#SBATCH --error=logs/checkm_per_sample_%j.err
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=24
#SBATCH --mem=32G

# === ENVIRONNEMENT ===
ENV_NAME="mahabio_env"
BASE_INPUT="/shared/home/asandri/MAHABIO_analysis/results/binning"
BASE_OUTPUT="/shared/home/asandri/MAHABIO_analysis/results/binning/checkm"
THREADS=8

# === ACTIVER CONDA ===
echo "🔁 Activation de l’environnement $ENV_NAME"
source ~/miniforge3/etc/profile.d/conda.sh
source activate "$ENV_NAME"

# === BINNER LISTE ===
BINNERS=("maxbin" "semibin" "vamb")

for BINNER in "${BINNERS[@]}"; do
    echo "🗂️  Traitement des échantillons dans $BINNER"
    SAMPLE_DIRS=($(find "$BASE_INPUT/$BINNER" -mindepth 1 -maxdepth 1 -type d))

    for SAMPLE_DIR in "${SAMPLE_DIRS[@]}"; do
        SAMPLE=$(basename "$SAMPLE_DIR")
        
        # Définir l'extension et le sous-répertoire selon le binner
        case "$BINNER" in
            maxbin)
                INDIR="$SAMPLE_DIR"
                EXT="fasta"
                ;;
            semibin)
                INDIR="$SAMPLE_DIR/output_bins"
                EXT="fa.gz"
                ;;
            vamb)
                INDIR="$SAMPLE_DIR/bins"
                EXT="fna"
                ;;
            *)
                echo "❌ Binner non reconnu : $BINNER"
                continue
                ;;
        esac

        OUTDIR="${BASE_OUTPUT}/${BINNER}/${SAMPLE}"
        if [ -d "$OUTDIR" ]; then
        echo "⚠️  Dossier de sortie déjà existant pour $BINNER / $SAMPLE. On saute."
        continue
    fi
    mkdir -p "$OUTDIR"

        echo "➡️  Lancement CheckM pour $BINNER / $SAMPLE"
        checkm lineage_wf -x "$EXT" --threads "$THREADS" "$INDIR" "$OUTDIR" > "$OUTDIR/checkm.log" 2>&1
        checkm qa --tab_table --out_format 2 "$OUTDIR/lineage.ms" "$OUTDIR" > "$OUTDIR/quality_summary.tsv"
        echo "✅ CheckM terminé pour $BINNER / $SAMPLE"
    done
done

echo "🎉 Tous les traitements CheckM par échantillon sont terminés."
