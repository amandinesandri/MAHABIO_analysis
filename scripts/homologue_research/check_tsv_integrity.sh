#!/bin/bash

# Usage : ./check_tsv_integrity.sh <file.tsv>

FILE="$1"

if [[ ! -f "$FILE" ]]; then
    echo "❌ Fichier introuvable : $FILE"
    exit 1
fi

echo "🔍 Vérification du nombre de colonnes dans : $FILE"
echo

# Nombre de colonnes sur chaque ligne
awk -F'\t' '{print NF}' "$FILE" | sort | uniq -c

# Test de cohérence
echo
EXPECTED=$(head -n 1 "$FILE" | awk -F'\t' '{print NF}')
ACTUAL=$(awk -F'\t' '{print NF}' "$FILE" | sort -nu | wc -l)

if [[ "$ACTUAL" -eq 1 ]]; then
    echo "✅ Fichier correct : $EXPECTED colonnes détectées uniformément."
    exit 0
else
    echo "⚠️ Attention : lignes avec nombre de colonnes variables détectées !"
    echo "👉 Lignes problématiques (≠ $EXPECTED colonnes) :"
    awk -F'\t' -v exp="$EXPECTED" 'NF != exp {print NR ": " $0}' "$FILE"
    exit 2
fi
