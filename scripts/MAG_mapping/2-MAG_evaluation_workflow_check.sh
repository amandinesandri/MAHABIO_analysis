#!/bin/bash
#SBATCH --job-name=MAG_evaluation_workflow_check
#SBATCH --output=logs/MAG_evaluation_workflow_check_%A.out
#SBATCH --error=logs/MAG_evaluation_workflow_check_%A.err
#SBATCH --time=72:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G

# Charger conda ou modules si besoin
# module load anaconda3
source activate mahabio_env

echo "=== Lancement du pipeline d'analyse MAHABIO ==="

# Lancer l'array pour dRep + GTDB-Tk + CompareM
sbatch --array=0-6 /shared/home/asandri/MAHABIO_analysis/scripts/MAG_mapping/MAG_evaluation.sh
array_job_id=$!

echo "Array lancé avec Job ID: ${array_job_id}"

# Attendre la fin de l'array
echo "=== En attente de la fin de l'array... ==="
scontrol show job ${array_job_id}
sacct --format=JobID,State

# Boucle pour vérifier régulièrement si l'array est terminé
while true; do
    state=$(sacct -j ${array_job_id} --format=State%20 -n | head -n1 | awk '{print $1}')
    
    if [[ "$state" == "COMPLETED" ]]; then
        echo "Tous les jobs sont terminés avec succès."
        break
    elif [[ "$state" == "FAILED" || "$state" == "CANCELLED" ]]; then
        echo "Erreur détectée dans l'array ! État: $state"
        exit 1
    else
        echo "Jobs encore en cours... état actuel: $state"
        sleep 300  # Attendre 5 minutes avant de re-checker
    fi
done

# Après la complétion de l'array, lancer la fusion
echo "=== Lancement du merge final des MAGs ==="
python3 /shared/home/asandri/MAHABIO_analysis/scripts/MAG_mapping/merge_MAG_info.py

echo "=== Pipeline complet ! Résultats fusionnés disponibles ==="
