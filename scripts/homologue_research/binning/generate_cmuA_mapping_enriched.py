# generate_cmuA_mapping_enriched.py

import os
import pandas as pd
import sys

def parse_header(header_line):
    """Parse a Prodigal/CheckM style header line."""
    parts = header_line.strip().split("#")
    if len(parts) < 2:
        return None, None, None, None, None

    coords = parts[1].strip().split()
    meta = parts[2].strip() if len(parts) > 2 else ""

    start, end, strand = (None, None, None)
    partial, gc_cont = (None, None)

    if len(coords) >= 3:
        start = coords[0]
        end = coords[1]
        strand = coords[2]

    meta_fields = dict()
    for item in meta.split(";"):
        if "=" in item:
            k, v = item.split("=")
            meta_fields[k.strip()] = v.strip()

    partial = meta_fields.get("partial", "NA")
    gc_cont = meta_fields.get("gc_cont", "NA")

    return start, end, strand, partial, gc_cont

def main(hits_tsv, faa_dir, output_file):
    hits = pd.read_csv(hits_tsv, sep="\t")

    records = []
    for idx, row in hits.iterrows():
        seq_id = f"seq_{idx+1}"
        binner = row['binner']
        sample = row['sample']
        bin_number = str(row['bin_number'])
        ref_bin = row['ref_bin']
        target_name = row['target_name']

        # Déterminer le chemin du fichier FAA
        faa_file = os.path.join(faa_dir, f"{ref_bin}_cmuA_hits.faa")

        start = end = strand = partial = gc_cont = "NA"

        if os.path.exists(faa_file):
            with open(faa_file) as f:
                for line in f:
                    if line.startswith(">") and target_name in line:
                        start, end, strand, partial, gc_cont = parse_header(line)
                        break

        records.append([seq_id, binner, sample, bin_number, ref_bin, target_name, start, end, strand, partial, gc_cont])

    enriched_df = pd.DataFrame(records, columns=[
        "seq_id", "binner", "sample", "bin_number", "ref_bin", "target_name",
        "start", "end", "strand", "partial", "gc_content"
    ])

    enriched_df.to_csv(output_file, sep="\t", index=False)
    print(f"✅ Mapping enrichi généré : {output_file}")

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python generate_cmuA_mapping_enriched.py hits_tsv faa_dir output_file")
        sys.exit(1)

    hits_tsv = sys.argv[1]
    faa_dir = sys.argv[2]
    output_file = sys.argv[3]

    main(hits_tsv, faa_dir, output_file)
