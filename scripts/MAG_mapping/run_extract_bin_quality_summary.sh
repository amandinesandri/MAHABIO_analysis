#!/bin/bash
#SBATCH --job-name=bin_quality_summary
#SBATCH --output=logs/bin_quality_summary.out
#SBATCH --error=logs/bin_quality_summary.err
#SBATCH --time=00:10:00
#SBATCH --mem=2G
#SBATCH --cpus-per-task=1

# Activation de l'environnement
source activate mahabio_env

# Lancement du script Python
python3 extract_bin_quality_summary.py
