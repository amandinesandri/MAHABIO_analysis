import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches

def plot_cmuA_positions_with_distances(tsv_file, output_file="cmuA_ORFs_plot.png"):
    df = pd.read_csv(tsv_file, sep="\t")
    df = df.sort_values(by="Start").reset_index(drop=True)

    fig, ax = plt.subplots(figsize=(12, 1 + len(df) * 0.4))

    y = 0
    for idx, row in df.iterrows():
        start = row["Start"]
        end = row["End"]
        strand = row["Strand"]
        length = row["Length"]
        orf_id = row["ORF_ID"]
        color = "skyblue" if strand == "+" else "salmon"
        direction = 1 if strand == "+" else -1
        width = end - start

        # Flèche représentant l'ORF
        arrow = patches.FancyArrow(
            start, y,
            width * direction, 0,
            width=0.3,
            head_length=min(500, width * 0.3),
            length_includes_head=True,
            color=color,
            edgecolor='black'
        )
        ax.add_patch(arrow)

        # Annotation du gène
        ax.text(start, y + 0.4, f"ORF {orf_id}\n{length} bp", fontsize=8)

        # Distance avec le précédent
        if idx > 0:
            prev_end = df.loc[idx - 1, "End"]
            distance = start - prev_end
            if distance > 0:
                mid = (prev_end + start) / 2
                ax.annotate(
                    f"{distance} bp",
                    xy=(mid, y),
                    xytext=(mid, y - 0.5),
                    ha='center',
                    fontsize=7,
                    arrowprops=dict(arrowstyle='->', lw=0.5, color='gray')
                )

        y += 1

    ax.set_xlim(0, df["End"].max() + 10000)
    ax.set_ylim(-1, y + 1)
    ax.set_xlabel("Position sur le contig (bp)")
    ax.set_yticks([])
    ax.set_title("Positions des CDS cmuA sur le contig Sj|ctg000000063")
    plt.tight_layout()
    plt.savefig(output_file)
    print(f"✅ Figure enregistrée : {output_file}")

if __name__ == "__main__":
    import sys
    if len(sys.argv) != 2:
        print("Usage : python plot_cmuA_orfs_with_distances.py cmuA_positions_enriched.tsv")
    else:
        plot_cmuA_positions_with_distances(sys.argv[1])
