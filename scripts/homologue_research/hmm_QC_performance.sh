#!/bin/bash
#SBATCH --job-name=QC_check_hmm
#SBATCH --output=logs/QC_check_hmm.out
#SBATCH --error=logs/QC_check_hmm.err
#SBATCH --time=01:00:00
#SBATCH --mem=4G

module load conda
source activate mahabio_env

cd /shared/home/asandri/MAHABIO_analysis
mkdir results/hmm/QC/
# Vérif consensus performance avec hmsearch sur fichier avec seq de ref mais je trouve pas ça si pertinent 

hmmsearch --tblout results/hmm/QC/QC_ref_hits.tbl results/hmm/cmuA.hmm results/hmm/cmuA_aligned.fasta
# générer une protéine représentative du modèle
hmmemit -c results/hmm/cmuA.hmm > results/hmm/QC/cmuA_consensus.faa

# générer plusieur séquences aléatoires représentatives du HMM
hmmemit -n 10 results/hmm/cmuA.hmm > results/hmm/QC/cmuA_10_random.faa

cd /shared/home/asandri/MAHABIO/data/
wget https://ftp.ncbi.nih.gov/refseq/release/bacteria/bacteria.1447.protein.faa.gz 

gunzip bacteria.1447.protein.faa.gz 

cat bacteria.1447.protein.faa> /shared/home/asandri/MAHABIO/data/refseq_bacteria.faa
cd /shared/home/asandri/MAHABIO_analysis


hmmsearch \
  --cpu 4 \
  --tblout results/hmm/QC/cmuA_refseq_ncbi.tbl \
  --domtblout results/hmm/QC/cmuA_refseq_ncbi.domtbl \
  results/hmm/cmuA.hmm \
  /shared/home/asandri/MAHABIO/data/refseq_bacteria.faa
