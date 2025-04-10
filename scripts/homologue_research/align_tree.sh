#!/bin/bash
#SBATCH --job-name=tree
#SBATCH --output=logs/tree_%j.out
#SBATCH --error=logs/tree_%j.err
#SBATCH --time=01:00:00
#SBATCH --mem=8G

source activate mahabio_env

mkdir ../../results/hmm/tree



#cat ../../results/hmm/hits/*_hits.faa /shared/home/asandri/MAHABIO/data/cmuaA_seq_prot.fasta > ../../results/hmm/tree/all_CmuA.fasta
#seqkit replace -p '\*' -r '' ../../results/hmm/tree/all_CmuA.fasta > ../../results/hmm/tree/all_CmuA_cleaned.faa
#mafft --auto ../../results/hmm/tree/all_CmuA_cleaned.faa > ../../results/hmm/tree/aligned_all_CmuA.fasta
#FastTree -lg ../../results/hmm/tree/aligned_all_CmuA.fasta > ../../results/hmm/tree/CmuA_tree.nwk

#mafft --auto /shared/home/asandri/MAHABIO/data/cmuaA_seq_prot.fasta > ../../results/hmm/tree/aligned_refseq_for_TEST_CmuA.fasta
#FastTree -lg ../../results/hmm/tree/aligned_refseq_for_TEST_CmuA.fasta > ../../results/hmm/tree/refseq_TEST_CmuA_tree.nwk



cat ../../results/hmm/hmm_matched/*_matched.faa  /shared/home/asandri/MAHABIO/data/cmuaA_seq_prot.fasta > ../../results/hmm/tree/all_CmuA.fasta
#seqkit replace -p '\*' -r '' ../../results/hmm/tree/all_CmuA.fasta > ../../results/hmm/tree/all_CmuA_cleaned.faa
mafft --auto ../../results/hmm/tree/all_CmuA.fasta > ../../results/hmm/tree/aligned_all_CmuA.fasta
FastTree -lg ../../results/hmm/tree/aligned_all_CmuA.fasta > ../../results/hmm/tree/CmuA_tree.nwk