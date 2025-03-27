import sys
import pandas as pd

def parse_domtblout(domtbl_file, evalue_cutoff=1e-5, coverage_cutoff=0.5):
    data = []
    with open(domtbl_file) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split()
            try:
                target_id = parts[0]
                query_id = parts[3]
                seq_len = int(parts[2])
                ali_start = int(parts[17])
                ali_end = int(parts[18])
                evalue = float(parts[12])
                score = float(parts[13])
                aln_len = abs(ali_end - ali_start) + 1
                coverage = aln_len / seq_len
                #if evalue <= evalue_cutoff and coverage >= coverage_cutoff:
                if True:
                    data.append([
                        target_id, query_id, evalue, score,
                        ali_start, ali_end, seq_len, aln_len, coverage
                    ])
            except (IndexError, ValueError):
                continue
    df = pd.DataFrame(data, columns=[
        "sequence_id", "domain_name", "evalue", "score",
        "align_start", "align_end", "sequence_length",
        "domain_length", "coverage"
    ])
    return df

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python parse_hmm_domtbl.py input.domtbl output.tsv")
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2]

    df = parse_domtblout(input_file)
    df.to_csv(output_file, sep="\t", index=False)
