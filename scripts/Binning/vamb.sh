#!/bin/bash
#SBATCH --job-name=vamb
#SBATCH --output=logs/vamb.out
#SBATCH --error=logs/vamb.err
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --array=0-2


# Activation de l'environnement conda contenant VAMB
source activate vamb_env  # ou le nom de ton env contenant VAMB

# Paramètres
contig_files=("/shared/home/asandri/MAHABIO/data/C_contigs_more_than_300bp.fasta" "/shared/home/asandri/MAHABIO/data/Sj_contigs_more_than_300bp.fasta" "/shared/home/asandri/MAHABIO/data/Fk_contigs_more_than_300bp.fasta")
BAM_DIR=("/shared/home/asandri/MAHABIO_analysis/results/coverage/C_contigs_more_than_300bp_sorted.bam" "/shared/home/asandri/MAHABIO_analysis/results/coverage/Sj_contigs_more_than_300bp_sorted.bam" "/shared/home/asandri/MAHABIO_analysis/results/coverage/Fk_contigs_more_than_300bp_sorted.bam")

OUTDIR="/shared/home/asandri/MAHABIO_analysis/results/binning/vamb"

mkdir -p ${OUTDIR}
cd ${OUTDIR}

for i in "${!contig_files[@]}"; do
     # Création du fichier de profil de couverture au format VAMB
     echo "Création du fichier .tsv de profils de couverture"
     vamb --outdir "$OUTDIR" --fasta "${contigs_files[$i]}" \
          --bamfiles "${BAM_DIR[$i]}" \
          --minfasta 200000 \
          --threads 16
done