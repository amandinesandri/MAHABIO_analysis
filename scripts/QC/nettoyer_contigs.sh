#!/bin/bash
#SBATCH --job-name=clean_contigs
#SBATCH --output=logs/clean_contigs%j.out
#SBATCH --error=logs/clean_contigs%j.err
#SBATCH --time=02:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=4


# Répertoire contenant les fichiers FASTA
DATA_DIR="/shared/home/asandri/MAHABIO/data"

# Motif des fichiers à traiter
FILE_PATTERN="*contigs_more_than_300bp.fasta"

# Caractères à supprimer ou remplacer dans les noms de contigs
# Ici, on supprime le caractère '|'
CLEAN_PATTERN='s/|//g'

# Boucle sur chaque fichier correspondant au motif
for fasta_file in "$DATA_DIR"/$FILE_PATTERN; do
    # Vérifie si le fichier existe
    if [[ -f "$fasta_file" ]]; then
        echo "Traitement du fichier : $fasta_file"

        # Détermine le nom du fichier de sortie
        cleaned_file="${fasta_file%.fasta}_cleaned.fasta"

        # Utilise awk pour traiter le fichier
        awk -v pattern="$CLEAN_PATTERN" '
            BEGIN { OFS = "" }
            /^>/ {
                # Nettoie le nom du contig
                header = substr($0, 2)
                gsub(/\|/, "", header)
                print ">" header
                next
            }
            {
                print
            }
        ' "$fasta_file" > "$cleaned_file"

        echo "Fichier nettoyé créé : $cleaned_file"
    else
        echo "Fichier introuvable : $fasta_file"
    fi
done
