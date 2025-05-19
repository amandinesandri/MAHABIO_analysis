import csv
import re
import sys

def parse_hits_table(hits_tsv):
    """Lit la table de hits pour récupérer les IDs ORF (ex: 1267)"""
    orf_ids = []
    with open(hits_tsv) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.strip().split()
            contig_orf = parts[4]  # ex: Sj|ctg000000063_1267
            match = re.search(r"_(\d+)$", contig_orf)
            if match:
                orf_ids.append(match.group(1))
    return orf_ids

def parse_gff(gff_file, orf_ids):
    """Extrait les infos de position des ORFs d’intérêt"""
    info = []
    with open(gff_file) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.strip().split("\t")
            if len(parts) < 9:
                continue
            contig = parts[0]
            start = parts[3]
            end = parts[4]
            strand = parts[6]
            attr = parts[8]
            match = re.search(r"ID=(?:[^_]+_)?(\d+)", attr)
            if match and match.group(1) in orf_ids:
                info.append([contig, match.group(1), start, end, strand])
    return info

def main(hits_table, gff_file, output_tsv):
    orf_ids = parse_hits_table(hits_table)
    annotations = parse_gff(gff_file, orf_ids)
    with open(output_tsv, "w", newline="") as out:
        writer = csv.writer(out, delimiter="\t")
        writer.writerow(["Contig", "ORF_ID", "Start", "End", "Strand"])
        writer.writerows(annotations)
    print(f"Infos de position écrites dans : {output_tsv}")

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage : python extract_orf_positions.py hits.tsv gff_file output.tsv")
        sys.exit(1)
    main(sys.argv[1], sys.argv[2], sys.argv[3])
