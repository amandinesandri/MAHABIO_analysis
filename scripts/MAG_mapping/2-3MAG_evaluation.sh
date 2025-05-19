#!/bin/bash
#SBATCH --job-name=mag_evaluation
#SBATCH --output=logs/mag_evaluation_%A_%a.out
#SBATCH --error=logs/mag_evaluation_%A_%a.err
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --array=0-6
source activate mahabio_env

# Liste des binners et samples (adapter si nécessaire)
binners=(maxbin semibin semibin semibin vamb vamb vamb )
samples=(C_contigs_more_than_300bp C_contigs_more_than_300bp Fk_contigs_more_than_300bp Sj_contigs_more_than_300bp C_contigs_more_than_300bp Fk_contigs_more_than_300bp Sj_contigs_more_than_300bp)

# Variables liées à l'array
binner=${binners[$SLURM_ARRAY_TASK_ID]}
sample=${samples[$SLURM_ARRAY_TASK_ID]}

# Définir les chemins d'entrée en fonction du binner
if [ "$binner" == "semibin" ]; then
    input_path="/shared/home/asandri/MAHABIO_analysis/results/binning/Affiliation/16S_detection/semibin/${sample}/*.fasta"
elif [ "$binner" == "maxbin" ]; then
    input_path="/shared/home/asandri/MAHABIO_analysis/results/binning/maxbin/${sample}/*.fasta"
elif [ "$binner" == "vamb" ]; then
    input_path="/shared/home/asandri/MAHABIO_analysis/results/binning/vamb/${sample}/bins/*.fna"
fi

# Définir les chemins de sortie
output_path="/shared/home/asandri/MAHABIO_analysis/results/binning/dRep/${binner}/${sample}"
gtdbtk_out="/shared/home/asandri/MAHABIO_analysis/results/binning/GTDBtk/${binner}/${sample}"
comparem_out="/shared/home/asandri/MAHABIO_analysis/results/binning/CompareM/${binner}/${sample}"

# # Vérifier si l'output dRep existe déjà
# if [ -d "$output_path" ]; then
#     echo "Analyse dRep déjà faite pour $binner/$sample --> SKIP dRep"
# else
#     Créer dossier de sortie
#     mkdir -p "$output_path"

#     # Charger l'environnement

#     echo "=== Lancement dRep pour $binner/$sample ==="
#     dRep dereplicate "$output_path" \
#         -g $input_path \
#         --completeness 50 \
#         --contamination 10 \
#         -pa 0.9 \
#         -sa 0.99 \
#         --processors 16
# fi



# # Maintenant lancer GTDB-Tk si pas déjà fait
# if [ ! -f "${gtdbtk_out}/gtdbtk.bac120.summary.tsv" ]; then
#     echo "=== Lancement GTDB-Tk pour $binner/$sample ==="
#     mkdir -p "$gtdbtk_out"
#     #source activate gtdbtk_env
#     gtdbtk classify_wf \
#         --genome_dir "${output_path}/dereplicated_genomes" \
#         --out_dir "$gtdbtk_out" \
#         --cpus 16
# else
#     echo "Classification GTDB-Tk déjà faite pour $binner/$sample --> SKIP"
# fi

#Lancer CompareM AAI si pas déjà fait
if [ ! -f "${comparem_out}/aai_results.tsv" ]; then
    echo "=== Lancement CompareM Amerged_MAG_infoAI pour $binner/$sample ==="
    mkdir -p "$comparem_out"
    #source activate comparem_env
    comparem aai_wf \
        "${output_path}/dereplicated_genomes" \
        "$comparem_out" \
        --cpus 16 \
        --file_ext ".fasta"
else
    echo "Analyse CompareM AAI déjà faite pour $binner/$sample --> SKIP"
fi

echo "=== FIN TACHE $binner/$sample ==="
