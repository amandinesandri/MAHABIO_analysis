#!/bin/bash
#SBATCH --job-name=bam_check
#SBATCH --output=logs/bam_check_%j.out
#SBATCH --error=logs/bam_check_%j.err
#SBATCH --time=08:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=1

module load samtools

# Fichier contenant les tailles de bins
MAG_SIZES_TSV="/shared/home/asandri/MAHABIO_analysis/results/MAG_mapping/MAG_sizes.tsv"

# Résultat
output="/shared/home/asandri/MAHABIO_analysis/results/MAG_mapping/bam_verification_report.tsv"
echo -e "fichier\tstatus\talignements\tcouverture_moyenne\tprofondeur_moyenne\ttaille_bin_bp\tnb_contigs" > "$output"

# Petite fonction pour récupérer la taille d'un bin à partir du TSV
get_bin_size() {
    local mag_name=$1
    grep -P "\t${mag_name}\t" "$MAG_SIZES_TSV" | awk -F'\t' '{print $4}'
}

# Fonction pour compter les contigs dans un fichier fasta correspondant
count_contigs() {
    local fasta_file=$1
    if [[ -f "$fasta_file" ]]; then
        if [[ "$fasta_file" == *.gz ]]; then
            gunzip -c "$fasta_file" | grep -c "^>"
        else
            grep -c "^>" "$fasta_file"
        fi
    else
        echo "0"
    fi
}

find /shared/home/asandri/MAHABIO_analysis/bam_files_full/ -name "*.bam" | sort | while read -r bam_file; do
    bam_basename=$(basename "$bam_file")
    fasta_candidate=$(echo "$bam_basename" | sed -E 's/_long\.bam$|_short\.bam$/.fna/') # Ajuste selon ton nommage
    fasta_path=$(find /shared/home/asandri/MAHABIO_analysis/results/binning/ -name "$fasta_candidate" 2>/dev/null | head -n 1)
    
    taille_bin="NA"
    nb_contigs="NA"
    
    if [[ -n "$fasta_path" ]]; then
        taille_bin=$(get_bin_size "$(basename "$fasta_path")")
        nb_contigs=$(count_contigs "$fasta_path")
    fi

    if samtools quickcheck "$bam_file"; then
        count=$(samtools view -c "$bam_file")
        if [ "$count" -gt 0 ]; then
            depth_file=$(mktemp)
            samtools depth "$bam_file" > "$depth_file"
            
            if [ -s "$depth_file" ]; then
                total_positions=$(awk 'END {print NR}' "$depth_file")
                total_depth=$(awk '{sum += $3} END {print sum}' "$depth_file")
                profondeur_moyenne=$(awk -v t="$total_positions" -v d="$total_depth" 'BEGIN {if (t>0) printf "%.2f", d/t; else print "0"}')
                couverture_moyenne=$total_positions
            else
                profondeur_moyenne=0
                couverture_moyenne=0
            fi

            echo -e "${bam_basename}\tOK\t${count}\t${couverture_moyenne}\t${profondeur_moyenne}\t${taille_bin}\t${nb_contigs}" >> "$output"
            rm "$depth_file"
        else
            echo -e "${bam_basename}\tOK_vide\t0\t0\t0\t${taille_bin}\t${nb_contigs}" >> "$output"
        fi
    else
        echo -e "${bam_basename}\tProblème\tNA\tNA\tNA\t${taille_bin}\t${nb_contigs}" >> "$output"
    fi
done
