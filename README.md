# README



Materials and Methods for MAHABIO Project

This document summarizes the bioinformatics analyses conducted in the MAHABIO project. Each section corresponds to a specific type of analysis with the tested parameters.

## 1.  Search data of referene and interest

### 1.1 Search Functional Gene  CmuA and HmtA homologues

* **Goal:** Detect homologous sequences of chloromethane degradation genes CmuA and HmtA in metagenomic assemblies and in reads before assembly. See adjacent genes. Check where are we located according phylogenetic tree 
* **Input:**Contigs as fasta files to translate into protein sequences using Prodigal, fasta file of reference CmuA protein sequences and hmtA gene + protein sequences
* **Tool/Package:** MAFFT, HMMER (hmmbuild, hmmsearch), transeq/seqkit/biopython, BLAST+, PSI-BLAST
* **Reference Sequences:** Provided protein FASTA files for CmuA, HmtA, and gene sequences for HmtA.
* **Steps:**
  * Align reference protein sequences (CmuA/HmtA) using MAFFT  
  * Build profile HMM using `hmmbuild`.
  * Predict protein sequences from contigs using seqkit
  * Search predicted proteins using `hmmsearch` and save the resulting HMM tables.
  * Extract regions that match to create fasta of them to analyse phylogenetic tree along with reference sequences used to build the profile.  
* **Parameters:**
  * MAFFT: `--auto`
  * HMMER: default parameters, `--tblout` and `domtblout` for summary output
* **Output:** fasta file of ref sequences alignement , HMM profile , List of detected homologous sequences with scores, identities, and e-values.


## Done / Next steps / To test / Follow-up

Done: 
- Check if consensus CmuaA is robust with 1) qulity control of hhmm research on initial ref seq and 2) research on NCBI results are consistant with the spcies of reference provided initially
- length of consensus: modèle cmuA.hmm contient 521 positions de correspondance (match states)

Immediate actions:
- hmmsearch sur quels critères pour identifier un hit?
- Identify UNQUE CmuA hits :  which criteria ?
- tables in excel Format for Françoise


Next steps:
- [ ] Homologous genz research for hmtA in contigs + raw data 
- [ ] Homologous gene research for CmuA in raw data 
- [ ] Calculate informative metrics such as ANI
- [ ] Make observations of which gene are surronding functionnal genes of interest  

- [ ] assign taxonomy to contigs where cmuA have been detected 
- [ ] Vérifier si les MAGs avec cmuA contiennent aussi des gènes de la publication
- [ ] Parsing HMM table to have a summary of relevant info about  detected hits. 
- [ ] Vérification de la présence du cluster cmuBCA
- [ ] Comparaison de la syntenie, présence sur plasmides vs chromosomes
- [ ] Alignement des hits avec un clade de référence (e.g. sensu stricto, cmuA-like, cmuA-anaérobie)
- [ ] Construction d’arbre phylogénétique pour classification fonctionnelle

On going:
- [ ] Very low quanitty of hits detected in Fk and Sj .. -> we keep them all ?
- [ ] No hit cmuA extrated as match from domtbl  so I remove coverage threshold 0.5 in extraction , 
- [ ] All samples CmuA hits sequences + CmuA reference sequences have been aigned and represented on a tree.



----------


### 1.2 Search the reference organism of chloromethane utilization: Methylobacterium Extorquence CM4
1.2 Search the reference organism of chloromethane utilization: Methylorubrum extorquens CM4
* **Goal:** Detect the presence of Methylorubrum extorquens CM4 (reference strain for chloromethane degradation via the cmu pathway) in metagenomic data, by identifying contigs, reads, or MAGs taxonomically assigned to this strain or very close taxa. Compare with the literature reference (e.g. Bringel et al., 2019) to check if we recovered this key strain in our dataset.
* **Input:**Assembled contigs in FASTA format, Translated protein sequences from contigs and long/short reads, MAGs if available, NCBI RefSeq genome/proteome for M. extorquens CM4 (downloaded from NCBI), Taxonomic database (e.g. GTDB, Kraken2 DB, or NCBI nt/nr)
* **Tool/Package:**  BLAST+ (blastn, blastp, tblastn), Kraken2 or Kaiju (for read-based or contig-level taxonomic assignment), GTDB-Tk (if working with MAGs), fastANI (optional: for genome similarity check), Prokka (optional: for MAG annotation)
* **Steps:**
  * Download the complete genome/proteome of M. extorquens CM4 from NCBI
  * Index the reference genome/proteome for BLAST
  * Search contigs or translated proteins using blastn or blastp against CM4 sequences
  * (Optional) Run Kraken2/Kaiju on reads or contigs for taxonomic assignment
  * Search MAGs with GTDB-Tk to see if one or more bins cluster near M. extorquens
  * Use ANI (fastANI) to compare any MAGs to CM4 genome for similarity >95%
* **Parameters:**
  * blastn: -evalue 1e-10 -outfmt 6
  * blastp: -evalue 1e-10 -outfmt 6
  * Kraken2: default database, or custom DB with CM4 included
  * GTDB-Tk: standard classification workflow on MAGs
  * fastANI: genome-to-genome identity threshold >95%
* **Output:** .
BLAST hit table of contigs or proteins matching CM4 sequences
List of contigs or MAGs assigned to M. extorquens or close relatives
ANI score (if available) of best MAG compared to CM4
Summary report: evidence of CM4 presence or absence in dataset



## 2. Metagenomic Binning (MAG Reconstruction)

### 2.1 Binning with MaxBin2

* **Goal:** Recover metagenome-assembled genomes (MAGs) from metagenomic assemblies.
* **Input:** Contigs (>300 bp), abundance/coverage files, reads (when available).
* **Tool/Package:** MaxBin2
* **Parameters:**
  * `-contig` input contig file
  * `-reads` or `-abund` (for coverage support)
  * `-min_contig_length 1000`
  * `-thread` for parallelization
* **Output:** Binned MAGs in FASTA format, bin statistics

## 3. 16S rRNA Gene Detection

### 3.1 Search for 16S sequences in MAGs/Contigs

* **Goal:** Identify and extract 16S rRNA sequences to determine taxonomic identity.
* **Input:** Contigs or MAGs
* **Tool/Package:** Barrnap or CheckM (16S marker analysis)
* **Parameters:** Default
* **Output:** 16S gene sequences with location and taxonomy (if available)

## 4. Assembly Quality Assessment

### 4.1 Quality Check with CheckM

* **Goal:** Assess completeness and contamination of MAGs
* **Input:** MAGs from MaxBin2
* **Tool/Package:** CheckM v1.0.18
* **Parameters:**
  * `lineage_wf` workflow
  * Input: directory with MAGs
  * Output: tabular summary of completeness and contamination
* **Output:** Table with completeness, contamination, marker lineage, etc.
