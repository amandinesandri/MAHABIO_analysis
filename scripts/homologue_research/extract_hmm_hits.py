import sys
from Bio import SeqIO

def extract_hits(tbl_file, faa_file, output_file, evalue_cutoff=1e-5):
    # Lire les identifiants de hits significatifs
    hits = set()
    with open(tbl_file) as f:
        for line in f:
            if line.startswith("#"):
                continue
            cols = line.strip().split()
            if len(cols) > 4 and float(cols[4]) <= evalue_cutoff:
                hits.add(cols[0])

    # Extraire les séquences correspondantes
    with open(output_file, "w") as out:
        for record in SeqIO.parse(faa_file, "fasta"):
            if record.id in hits:
                SeqIO.write(record, out, "fasta")

if __name__ == "__main__":
    extract_hits(sys.argv[1], sys.argv[2], sys.argv[3])
