import os
import pandas as pd
import sys

def parse_hmm_tbl(tbl_file):
    """Parse a hmmsearch .tbl file and extract hits."""
    hits = []
    with open(tbl_file) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.strip().split()
            try:
                target = parts[0]
                query = parts[2]
                evalue = float(parts[4])
                score = float(parts[5])
                hits.append((target, query, evalue, score))
            except (IndexError, ValueError):
                continue
    return hits

def main(input_dir, output_file, evalue_cutoff=None, score_cutoff=None):
    records = []
    for root, dirs, files in os.walk(input_dir):
        for file in files:
            if file.endswith("_cmuA.tbl"):
                tbl_path = os.path.join(root, file)
                rel_path = os.path.relpath(tbl_path, input_dir)
                parts = rel_path.split(os.sep)

                if len(parts) >= 3:
                    binner = parts[0]
                    sample = parts[1]
                    ref_bin = file.replace("_cmuA.tbl", "")
                else:
                    binner = "unknown"
                    sample = "unknown"
                    ref_bin = file.replace("_cmuA.tbl", "")

                # Extraire bin_number proprement depuis ref_bin
                if "_" in ref_bin:
                    bin_number = ref_bin.split("_")[-1]
                else:
                    bin_number = ref_bin  # fallback

                hits = parse_hmm_tbl(tbl_path)
                for target, query, evalue, score in hits:
                    if (evalue_cutoff is None or evalue <= evalue_cutoff) and \
                       (score_cutoff is None or score >= score_cutoff):
                        records.append([binner, sample, bin_number, ref_bin, target, query, evalue, score])

    df = pd.DataFrame(records, columns=["binner", "sample", "bin_number", "ref_bin", "target_name", "query_name", "evalue", "score"])
    df.to_csv(output_file, sep="\t", index=False)
    print(f"✅ Résultats sauvegardés dans {output_file} ({len(df)} hits trouvés)")

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python parse_cmuA_hits.py input_dir output_file [evalue_cutoff] [score_cutoff]")
        sys.exit(1)

    input_dir = sys.argv[1]
    output_file = sys.argv[2]

    evalue_cutoff = float(sys.argv[3]) if len(sys.argv) > 3 else None
    score_cutoff = float(sys.argv[4]) if len(sys.argv) > 4 else None

    main(input_dir, output_file, evalue_cutoff, score_cutoff)