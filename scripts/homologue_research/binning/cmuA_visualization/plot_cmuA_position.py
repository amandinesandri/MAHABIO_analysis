# plot_cmuA_positions.py
import matplotlib.pyplot as plt
import re

def parse_tbl(tbl_path):
    hits = []
    with open(tbl_path) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.strip().split()
            contig = fields[0]
            evalue = float(fields[4])
            score = float(fields[5])
            desc = " ".join(fields[18:])
            match = re.search(r"# (\d+) # (\d+) # ([\-1]+)", desc)
            if match:
                start, end, strand = int(match[1]), int(match[2]), match[3]
                hits.append((contig, start, end, strand, score, evalue))
    return hits

def plot_hits(hits, output_file):
    fig, ax = plt.subplots(figsize=(10, len(hits)))
    y = 0
    for contig, start, end, strand, score, evalue in hits:
        width = abs(end - start)
        color = 'skyblue' if strand == "1" else 'salmon'
        ax.barh(y, width, left=min(start, end), height=0.6, color=color, edgecolor='black')
        ax.text(min(start, end), y + 0.2, f"{contig}\nScore: {score:.1f}", fontsize=8)
        y += 1
    ax.set_xlabel("Position sur le contig")
    ax.set_yticks([])
    ax.set_title("Occurrences du gène cmuA sur bin")
    plt.tight_layout()
    plt.savefig(output_file)
    print(f"Figure saved to {output_file}")

# Exemple d’utilisation
if __name__ == "__main__":
    import sys
    tbl_file = sys.argv[1]
    output_file = sys.argv[2]
    hits = parse_tbl(tbl_file)
    plot_hits(hits, output_file)
