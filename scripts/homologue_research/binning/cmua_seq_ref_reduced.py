from Bio import SeqIO

# === Paramètres ===
input_fasta = "/shared/home/asandri/MAHABIO/data/cmuaA_seq_prot.fasta"  # Remplace par ton chemin exact si besoin
output_fasta = "/shared/home/asandri/MAHABIO/data/cmuaA_seq_prot_filtered.fasta"
keywords = ["Hyphomicrobium", "corrinoid", "methyltransferase", "Rhodobacteraceae", "Methylobacterium"]

# === Lecture et filtrage ===
filtered_seqs = []
for record in SeqIO.parse(input_fasta, "fasta"):
    header = record.description
    if any(kw.lower() in header.lower() for kw in keywords):
        # Nettoyage du header pour l'output
        clean_id = header.replace(" ", "_").replace(",", "").replace(";", "")
        clean_id = clean_id.split()[0]  # On garde le 1er mot comme identifiant
        record.id = f"ref_{clean_id}"
        record.description = ""  # Supprime la description pour garder juste l'ID
        filtered_seqs.append(record)

# === Écriture du fichier FASTA filtré ===
SeqIO.write(filtered_seqs, output_fasta, "fasta")
print(f"✅ {len(filtered_seqs)} séquences sauvegardées dans : {output_fasta}")
