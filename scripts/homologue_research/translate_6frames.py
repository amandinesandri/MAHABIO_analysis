# scripts/translate_6frames.py
import sys
from Bio import SeqIO
from Bio.Seq import Seq

def translate_6frames(input_fasta, output_fasta):
    with open(output_fasta, "w") as out:
        for record in SeqIO.parse(input_fasta, "fasta"):
            seq = record.seq
            for strand, nuc in [(+1, seq), (-1, seq.reverse_complement())]:
                for frame in range(3):
                    trans = nuc[frame:].translate(to_stop=False)
                    header = f"{record.id}_strand{strand}_frame{frame+1}"
                    out.write(f">{header}\n{str(trans)}\n")

if __name__ == "__main__":
    input_fasta = sys.argv[1]
    output_fasta = sys.argv[2]
    translate_6frames(input_fasta, output_fasta)
