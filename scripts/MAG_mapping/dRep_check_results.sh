#!/bin/bash
#SBATCH --job-name=drep_check
#SBATCH --output=logs/drep_check.out
#SBATCH --error=logs/drep_check.err
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:30:00

source activate mahabio_env
# Dossier de base
BASE_DIR=~/MAHABIO_analysis/results/binning/dRep

# Script Python inline
PY_SCRIPT=$(cat << 'EOF'
import pandas as pd
import sys
from pathlib import Path

base = Path(sys.argv[1])
report = []

for path in base.rglob('Chdb.csv'):
    drep_dir = path.parent.parent
    name = drep_dir.relative_to(base)
    
    chdb = pd.read_csv(path)
    filtered = chdb[(chdb['Completeness'] >= 50) & (chdb['Contamination'] <= 10)]
    
    derep_file = drep_dir / 'dereplicated_genomes'
    if derep_file.exists():
        with open(derep_file) as f:
            derep_genomes = set(line.strip() for line in f)
        derep_filtered = filtered[filtered['genome'].isin(derep_genomes)]
    else:
        derep_filtered = pd.DataFrame()
    
    out_file = drep_dir / f'{name.name}_filtered_genomes.tsv'
    derep_filtered.to_csv(out_file, sep='\t', index=False)
    
    report.append({
        'sample_binner': str(name),
        'dereplicated': len(derep_genomes) if derep_file.exists() else 0,
        'derep_passed_filters': len(derep_filtered)
    })

summary = pd.DataFrame(report)
summary.to_csv(base / 'summary_dRep_filters.tsv', sep='\t', index=False)
EOF
)

# Lancer le script Python
python3 -c "$PY_SCRIPT" "$BASE_DIR"
