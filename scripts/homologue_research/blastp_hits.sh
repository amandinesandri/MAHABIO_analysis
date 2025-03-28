#!/bin/bash
#SBATCH --job-name=blastp_hits
#SBATCH --output=logs/blastp_%x_%j.out
#SBATCH --error=logs/blastp_%x_%j.err
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G

# ⚙️ Chargement de l'environnement (modifie selon ton cas)
source activate mahabio_env  # Assure-toi que blastp est dedans

# 📁 Répertoires et fichiers
INPUT_FASTA=$1  
BASENAME=$(basename "$INPUT_FASTA" .faa)
OUTDIR=../../results/hmm/blast_hits
mkdir -p "$OUTDIR"

# 🚀 Lancement du BLASTp (ici version en ligne avec base 'nr')
blastp \
  -query "$INPUT_FASTA" \
  -db nr \
  -remote \
  -out "$OUTDIR/${BASENAME}_blast.tsv" \
  -evalue 1e-5 \
  -num_threads 4 \
  -max_target_seqs 1 \
  -outfmt '6 qseqid sseqid pident evalue bitscore staxids sscinames sskingdoms'

echo "BLAST terminé pour $INPUT_FASTA"
