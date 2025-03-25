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
  * Predict protein sequences from contigs using transeq/seqkit/biopython
  * Search predicted proteins using `hmmsearch` or PSI-BLAST with custom PSSM.
* **Parameters:**
  * MAFFT: `--auto`
  * HMMER: default parameters, `--tblout` for summary output
  * PSI-BLAST: `-num_iterations 3`, `-evalue 1e-5`, `-outfmt 6`
* **Output:** fasta file of ref sequences alignement , HMM profile, List of detected homologous sequences with scores, identities, and e-values.


### 1.2 Search the reference organism of chloromethane utilization: Methylobacterium Extorquence CM4

* **Goal:** 
* **Input:**
* **Tool/Package:** 
* **Reference Sequences:** 
* **Steps:**

* **Parameters:**

* **Output:** .




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

## General Metadata for Each Analysis

* **Dataset Used:** (e.g., `C_contigs_more_than_300bp.fasta`, `Fk_proteins.faa`, etc.)
* **Date of Execution:** Indicate the date the analysis was run.
* **Software Version:** Include version numbers (e.g., MAFFT 7.505, HMMER 3.4)
* **Environment:** Local (Ubuntu 20.04), Cluster (SLURM-based), KBase (cloud platform)

> This document will be updated iteratively as new analyses are performed or refined.

Feel free to suggest additions or reorganizations depending on future developments.




