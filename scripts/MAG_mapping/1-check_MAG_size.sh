#!/bin/bash
#SBATCH --job-name=check_MAG_size_filtered
#SBATCH --output=logs/check_MAG_size_filtered_%j.out
#SBATCH --error=logs/check_MAG_size_filtered_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:30:00

DETAILS_FILE="/shared/home/asandri/MAHABIO_analysis/results/MAG_mapping/MAG_sizes.tsv"
SUMMARY_FILE="/shared/home/asandri/MAHABIO_analysis/results/MAG_mapping/bins_counts_per_sample-binner.csv"
echo -e "Sample;Binner;N_bins;Total_size_Mbp;Mean_size_Mbp;Median_size_Mbp;Max_size_Mbp;Min_size_Mbp" > "$SUMMARY_FILE"

echo -e "ref_bin\tbinner\tSample\tbin_number\tfilename\tSize_bp\tSize_Mbp" > "$DETAILS_FILE"

declare -A sizes_map

# Fonction pour extraire le bin_number et écrire la ligne détaillée
add_entry () {
    local binner=$1
    local sample=$2
    local filepath=$3

    filename=$(basename "$filepath")
    # Supprimer les extensions classiques .fasta, .fa, .fa.gz, .fna, etc.
    if [[ "$filename" == *.fa.gz ]]; then
        name_without_ext="${filename%.fa.gz}"
    elif [[ "$filename" == *.fasta ]]; then
        name_without_ext="${filename%.fasta}"
    elif [[ "$filename" == *.fa ]]; then
        name_without_ext="${filename%.fa}"
    elif [[ "$filename" == *.fna ]]; then
        name_without_ext="${filename%.fna}"
    else
        name_without_ext="$filename"  # Pas d’extension reconnue
    fi

    if [[ "$binner" == "maxbin" ]]; then
        # Pour MaxBin : extraction du numéro spécifique
        bin_number=$(echo "$filename" | sed -E 's/.*maxbin\.0*([0-9]+)\.fasta/\1/')
    else
        # Pour VAMB et SemiBin : extraction générique
        bin_number=$(echo "$filename" | sed -E 's/[^0-9]*([0-9]+).*/\1/')
    fi

    # Taille brute en bp
    if [[ "$filepath" == *.gz ]]; then
        size_bp=$(gunzip -c "$filepath" | grep -v "^>" | tr -d '\n' | wc -c)
    else
        size_bp=$(grep -v "^>" "$filepath" | tr -d '\n' | wc -c)
    fi

    size_mbp=$(awk "BEGIN {printf \"%.2f\", $size_bp / 1000000}")
    ref_bin="${sample}_${binner}_${bin_number}"

    echo -e "${ref_bin}\t${binner}\t${sample}\t${bin_number}\t${name_without_ext}\t${size_bp}\t${size_mbp}" >> "$DETAILS_FILE"

    key="${sample}_${binner}"
    sizes_map["$key"]+="${size_mbp} "
}


# === VAMB ===
for sample in C_contigs_more_than_300bp Fk_contigs_more_than_300bp Sj_contigs_more_than_300bp; do
    for fna in /shared/home/asandri/MAHABIO_analysis/results/binning/vamb/${sample}/bins/*.fna; do
        [[ -e "$fna" ]] && add_entry "vamb" "$sample" "$fna"
    done
done

# === MAXBIN ===
for sample in C_contigs_more_than_300bp; do
    for fasta in /shared/home/asandri/MAHABIO_analysis/results/binning/maxbin/${sample}/*.fasta; do
        [[ -e "$fasta" ]] && add_entry "maxbin" "$sample" "$fasta"
    done
done

# === SEMIBIN ===
for sample in C_contigs_more_than_300bp Fk_contigs_more_than_300bp Sj_contigs_more_than_300bp; do
    for fa in /shared/home/asandri/MAHABIO_analysis/results/binning/semibin/${sample}/output_bins/*.fa.gz; do
        [[ -e "$fa" ]] && add_entry "semibin" "$sample" "$fa"
    done
done

# === Statistiques résumé ===
for key in "${!sizes_map[@]}"; do
    IFS="_" read sample binner <<< "$key"
    sizes=(${sizes_map[$key]})
    count=${#sizes[@]}
    total=0
    for val in "${sizes[@]}"; do
        total=$(awk -v t=$total -v v=$val 'BEGIN {print t + v}')
    done
    mean=$(awk -v t=$total -v c=$count 'BEGIN {printf "%.2f", t / c}')
    sorted=($(printf '%s\n' "${sizes[@]}" | sort -n))
    if (( count % 2 == 1 )); then
        median=${sorted[$((count / 2))]}
    else
        m1=${sorted[$((count / 2 - 1))]}
        m2=${sorted[$((count / 2))]}
        median=$(awk -v a=$m1 -v b=$m2 'BEGIN {printf "%.2f", (a + b) / 2}')
    fi
    max=${sorted[-1]}
    min=${sorted[0]}
    echo "$sample;$binner;$count;$total;$mean;$median;$max;$min" >> "$SUMMARY_FILE"
done

echo "✅ Résultats écrits dans :"
echo "- $DETAILS_FILE"
echo "- $SUMMARY_FILE"
