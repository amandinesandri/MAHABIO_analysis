import sys
from Bio import SeqIO

def parse_domtbl(domtbl_file, evalue_cutoff=1e-5, coverage_min=0.5):
    hits = {}
    with open(domtbl_file) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split()
            target_id = parts[0]
            query_id = parts[3]
            full_seq_len = int(parts[2])
            ali_start = int(parts[17])
            ali_end = int(parts[18])
            evalue = float(parts[12])
            aln_len = abs(ali_end - ali_start) + 1
            coverage = aln_len / full_seq_len

            #if evalue <= evalue_cutoff and coverage >= coverage_min:
            if True:
                hits[target_id] = (ali_start, ali_end)
    return hits

def extract_regions(fasta_file, hits_dict, output_file):
    with open(output_file, "w") as out:
        for record in SeqIO.parse(fasta_file, "fasta"):
            if record.id in hits_dict:
                start, end = hits_dict[record.id]
                if start < end:
                    sub = record.seq[start-1:end]
                else:
                    sub = record.seq[end-1:start]
                record.seq = sub
                SeqIO.write(record, out, "fasta")

if __name__ == "__main__":
    domtbl = sys.argv[1]
    fasta = sys.argv[2]
    output = sys.argv[3]
    hits = parse_domtbl(domtbl)
    extract_regions(fasta, hits, output)
