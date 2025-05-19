# tbl_to_bed.py
import re

def tbl_to_bed(tbl_file, bed_file):
    with open(tbl_file) as fin, open(bed_file, "w") as fout:
        for line in fin:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.strip().split()
            contig = fields[0]
            desc = " ".join(fields[18:])
            match = re.search(r"# (\d+) # (\d+) # ([\-1]+)", desc)
            if match:
                start = int(match[1])
                end = int(match[2])
                strand = '+' if match[3] == '1' else '-'
                fout.write(f"{contig}\t{start}\t{end}\tcmuA\t0\t{strand}\n")

if __name__ == "__main__":
    import sys
    tbl_to_bed(sys.argv[1], sys.argv[2])
