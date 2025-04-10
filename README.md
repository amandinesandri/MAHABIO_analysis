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
