#!/bin/bash
#SBATCH --job-name=align_full_array
#SBATCH --output=logs/align_full_array/align_full_array_%A_%a.out
#SBATCH --error=logs/align_full_array/align_full_array_%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=24:00:00
#SBATCH --array=0-389

set -x

echo "✅ Job SLURM démarré : Task ID $SLURM_ARRAY_TASK_ID"

# Modules nécessaires
echo "🧬 Activation de l'environnement Conda"
source activate mahabio_env

echo "🔧 Chargement des modules minimap2, bwa, samtools"
module load minimap2
module load bwa
module load samtools

# Variables
MAG_DIR="/shared/home/asandri/MAHABIO_analysis/results/binning/"
READ_DIR="/shared/home/asandri/MAHABIO/data/"
BAM_DIR="/shared/home/asandri/MAHABIO_analysis/bam_files_full/"
mkdir -p $BAM_DIR

echo "📂 Construction de la liste des MAGs"
MAGS=(
$(ls ${MAG_DIR}/maxbin/C_contigs_more_than_300bp/*.fasta)
$(ls ${MAG_DIR}/semibin/*/output_bins/*.fa.gz)
$(ls ${MAG_DIR}/vamb/*/bins/*.fna)
)

MAG=${MAGS[$SLURM_ARRAY_TASK_ID]}
MAG_NAME=$(basename "$MAG")
echo "🔍 MAG sélectionné : $MAG_NAME"

# Détection méthode/échantillon
if [[ $MAG == *"maxbin"* ]]; then
    method="maxbin"
    sample="C"
elif [[ $MAG == *"semibin"* ]]; then
    method="semibin"
    [[ $MAG == *"C_contigs"* ]] && sample="C"
    [[ $MAG == *"Fk_contigs"* ]] && sample="Fk"
    [[ $MAG == *"Sj_contigs"* ]] && sample="Sj"
elif [[ $MAG == *"vamb"* ]]; then
    method="vamb"
    [[ $MAG == *"C_contigs"* ]] && sample="C"
    [[ $MAG == *"Fk_contigs"* ]] && sample="Fk"
    [[ $MAG == *"Sj_contigs"* ]] && sample="Sj"
fi

echo "📌 Méthode : $method | Échantillon : $sample"

# Préparation reads
short_R1="${READ_DIR}/${sample}_R1.fastq.gz"
short_R2="${READ_DIR}/${sample}_R2.fastq.gz"
long_reads="${READ_DIR}/${sample}.fastq.gz"

echo "📥 Reads short : $short_R1 / $short_R2"
echo "📥 Reads long  : $long_reads"

# Création des dossiers
mkdir -p ${BAM_DIR}/${method}/${sample}/long
mkdir -p ${BAM_DIR}/${method}/${sample}/short

# Décompression temporaire si besoin
echo "📦 Détection du format du MAG"
TEMP_FASTA=""
if [[ $MAG == *.fa.gz ]]; then
    echo "🗜️ Décompression de $MAG"
    gunzip -c "$MAG" > "temp_${MAG_NAME%.fa.gz}.fa"
    TEMP_FASTA="temp_${MAG_NAME%.fa.gz}.fa"
elif [[ $MAG == *.fasta || $MAG == *.fna ]]; then
    TEMP_FASTA="$MAG"
else
    echo "❌ Extension non reconnue pour $MAG — SKIP"
    exit 1
fi

BASENAME=$(basename "$TEMP_FASTA" .fa)
BASENAME=$(basename "$BASENAME" .fasta)
BASENAME=$(basename "$BASENAME" .fna)

# Mapping long reads
echo "🔄 Mapping des long reads avec minimap2"
minimap2 -ax map-ont "$TEMP_FASTA" "$long_reads" | samtools sort -@8 -o "${BAM_DIR}/${method}/${sample}/long/${BASENAME}_long.bam"
samtools index "${BAM_DIR}/${method}/${sample}/long/${BASENAME}_long.bam"

# Mapping short reads
echo "🔄 Indexation BWA"
bwa index "$TEMP_FASTA"
echo "🔄 Mapping des short reads avec BWA MEM"
bwa mem -t 8 "$TEMP_FASTA" "$short_R1" "$short_R2" | samtools sort -@8 -o "${BAM_DIR}/${method}/${sample}/short/${BASENAME}_short.bam"
samtools index "${BAM_DIR}/${method}/${sample}/short/${BASENAME}_short.bam"

# Nettoyage
if [[ $TEMP_FASTA == temp_* ]]; then
    echo "🧹 Nettoyage : suppression $TEMP_FASTA"
    rm "$TEMP_FASTA"
fi

echo "✅ Terminé : MAG $MAG_NAME aligné et BAMs créés pour $sample (TaskID=$SLURM_ARRAY_TASK_ID)"
