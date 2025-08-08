# Multi-omics Analysis Results for *Purpureocillium lilacinum* PLA-C1

This repository contains the outputs and intermediate results from the **Applied Genomics** analysis of *Purpureocillium lilacinum* PLA-C1, covering genome quality assessment, functional annotation, orthology inference, phylogenetic reconstruction, and transcriptomic differential expression.

---

## 00 – Input Data  
**Folder:** `00_Input_data/`  
Contains the raw genomic and transcriptomic input files used across analyses.  

- `genomic.gff` – Genome feature annotations.  
  *Used in:* transcript extraction (gffread), CAZyme annotation.  

- `GCF_023168085.1_PurlilCBS_1.0_genomic.zip` – Genomic assembly in FASTA format.  
  *Used in:* BUSCO genome completeness assessment.  

- `protein.faa` – Predicted protein sequences.  
  *Used in:* CAZyme annotation, OrthoVenn3 orthology analysis.  

- `rna.fna` – Transcript sequences.  
  *Used in:* transcriptome reference for RNA-Seq quantification.  

---

## 01 – BUSCO
**Folder:** `01_BUSCO/`  
Assessment of genome completeness using **BUSCO** (Benchmarking Universal Single-Copy Orthologs).

- `Galaxy14-[Busco on data 7_short summary].txt` – Summary table with completeness scores.

**Tool & Parameters:**  
BUSCO in `genome` mode with **auto-lineage** selection.

---

## 02 – CAZy Annotation
**Folder:** `02_CAZy_annotation/`  
Prediction of Carbohydrate-Active Enzymes using **dbCAN**.

- `CAZyme.pep` – Protein sequences identified as CAZymes.
- `h.out`, `hmmer.out` – HMMER search outputs against CAZy HMM database.
- `overview.txt` – Summary of CAZyme predictions.

**Tool & Parameters:**  
dbCAN3 with HMMER search, e-value < 1e-15, coverage > 0.35.

---

## 03 – Orthology Analysis
**Folder:** `03_Orthology_analysis/`  
Inference of orthologous groups with **OrthoFinder**.

- `Co-orthologs.txt` – Co-ortholog relationships.
- `Inparalogs.txt` – In-paralog gene lists.
- `Orthogroups.txt` – Orthogroup assignments.
- `Orthologs.txt` – Orthologous gene pairs.
- `Multiple_sequence_alignment.fasta` – Alignment of orthologous sequences.
  
- **figures/**
  - `jVenn_chart.png` – overlap of orthogroups between the three species.
  - `UpSetJS.png` – umber of shared and unique orthogroups among species.
  - `result_statics.png` – overlaps, total clusters, single-copy clusters, protein counts, and singleton percentage.
  - `table_data.png` – number of proteins, clusters, and singletons per species.
  - `pairwise_heatmap.png` – visual similarity in cluster content between species (darker red = higher similarity).

**Pipeline:**
1. OrthoFinder (orthogroup inference)  
2. MAFFT (alignment)  
3. Gblocks (conserved block selection)  
4. AMAS (concatenation)  

---

## 04 – Phylogeny
**Folder:** `04_Phylogeny/`  
Reconstruction of species relationships using single-copy orthologs.

- `f1a8dba68a6f43b189057b437429746f-fasta-tree.png` – Final phylogenetic tree image.  

**Pipeline:**
1. OrthoFinder → Single-copy orthologues.
2. MAFFT → Multiple sequence alignment.
3. Gblocks → Conserved block selection.
4. AMAS → Concatenated supermatrix.
5. Tree inference with MEGA, FastTree, RAxML, or MrBayes.

---

## 05 – Transcriptome
**Folder:** `05_Transcriptome/`  
Transcript abundance estimation and differential expression.

- `Galaxy11-[gffread on data 7 and data 10_pep.fa].fasta` – Transcript sequences extracted from GFF.
- *(to be added)* `rna_counts.tsv` – Transcript counts from Salmon quant.
- *(to be added)* `deseq2_results.csv` – Differentially expressed genes from DESeq2.

**Pipeline:**
1. **Salmon quant** – Mapping RNA-Seq reads to transcripts.
2. **tximport** – Generate `rna_counts.tsv`.
3. **DESeq2** – Identify differentially expressed genes.

---
