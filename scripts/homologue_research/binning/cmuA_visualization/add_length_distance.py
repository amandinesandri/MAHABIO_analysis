import pandas as pd
import sys

def compute_lengths_and_distances(input_tsv, output_tsv):
    df = pd.read_csv(input_tsv, sep="\t")

    # Convertir en entiers
    df["Start"] = df["Start"].astype(int)
    df["End"] = df["End"].astype(int)

    # Calculer la longueur
    df["Length"] = df["End"] - df["Start"] + 1

    # Trier par contig et position de début
    df = df.sort_values(by=["Contig", "Start"]).reset_index(drop=True)

    # Calculer les distances au précédent sur le même contig
    distances = []
    prev_pos = {}
    for idx, row in df.iterrows():
        contig = row["Contig"]
        start = row["Start"]
        if contig in prev_pos:
            distances.append(start - prev_pos[contig])
        else:
            distances.append(None)  # pas de précédent
        prev_pos[contig] = start
    df["Distance_to_previous"] = distances

    # Sauvegarde
    df.to_csv(output_tsv, sep="\t", index=False)
    print(f"✅ Fichier enrichi : {output_tsv}")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage : python add_length_and_distance.py input.tsv output.tsv")
        sys.exit(1)
    compute_lengths_and_distances(sys.argv[1], sys.argv[2])
