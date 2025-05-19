#!/bin/bash
#SBATCH --job-name=check_MAG_size_stats
#SBATCH --output=logs/check_MAG_size_stats_%j.out
#SBATCH --error=logs/check_MAG_size_stats_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:30:00

source activate mahabio_env

python 2-1-MAG_size_stat.py 