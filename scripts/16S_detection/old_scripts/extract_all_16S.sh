#!/bin/bash
#SBATCH --job-name=prep_extract_16S
#SBATCH --output=logs/extract_all_16S_%j.out
#SBATCH --error=logs/extract_all_16S_%j.err
#SBATCH --time=00:20:00
#SBATCH --mem=2G
#SBATCH --cpus-per-task=1

# ------------------ ⏱️ DÉBUT ------------------
START_TIME=$(date +%s)
echo "⏱️ Job lancé le : $(date '+%Y-%m-%d %H:%M:%S')"
# ------------------------------------------------

echo "📦 Préparation des extractions 16S..."

ROOT="/shared/home/asandri/MAHABIO_analysis/results/binning/checkm"
SAMPLES_LIST="samples_list.txt"

# Répertoire d'accueil des logs
mkdir -p logs

# Générer automatiquement la liste
> $SAMPLES_LIST
for binner in semibin maxbin vamb; do
    if [ -d "$ROOT/$binner" ]; then
        for sample in $(ls "$ROOT/$binner"); do
            echo "$binner/$sample" >> $SAMPLES_LIST
        done
    fi
done

# Combien d'échantillons ?
N_SAMPLES=$(wc -l < $SAMPLES_LIST)
N_SAMPLES=$((N_SAMPLES - 1)) # car --array commence à 0

echo "📝 Liste des samples :"
cat $SAMPLES_LIST
echo "📈 Nombre d'échantillons trouvés : $((N_SAMPLES+1))"
echo "🚀 Lancement du job array..."

# Lancer l'array
sbatch --array=0-${N_SAMPLES} extract_16s.sh $SAMPLES_LIST

echo "✅ Tout est lancé proprement."

# ------------------ ✅ FIN ------------------
END_TIME=$(date +%s)
DURATION=$((END_TIME - START_TIME))
echo "✅ Job terminé le : $(date '+%Y-%m-%d %H:%M:%S')"

hours=$((DURATION / 3600))
minutes=$(((DURATION % 3600) / 60))
seconds=$((DURATION % 60))
echo "⏳ Durée totale : ${hours}h ${minutes}min ${seconds}s"
# ------------------------------------------------
