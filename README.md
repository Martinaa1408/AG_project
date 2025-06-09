# De novo Genome Assembly and Variant Analysis of Cupriavidus necator H16

## Overview

This project focuses on the **de novo assembly and variant analysis** of an industrial strain of *Cupriavidus necator* H16, a model bacterium for polyhydroxyalkanoate (PHA) bioplastic production. The aim is to generate a high-quality genome assembly and identify genetic variants (SNPs, indels, CNVs) compared to the reference strain (DSM 428), with a focus on PHA biosynthetic genes (phaA, phaB, phaC).

---

## Objectives

-Generate a de novo genome assembly of an industrial strain of *C. necator* H16.  
-Compare this assembly with the reference genome to identify genomic variants.  
-Annotate the variants and assess their potential impact on PHA production.

---

## Materials and Methods

### 1️⃣ Sample Collection and DNA Extraction  
- Bacterial strain: industrial strain of *C. necator* H16.  
- DNA extraction using commercial kits.  
- DNA quality check: Nanodrop (A260/280), Qubit, agarose gel electrophoresis.

### 2️⃣ Sequencing  
- **Long-read sequencing**: PacBio HiFi reads (~100X coverage).  
- **Short-read sequencing**: Illumina paired-end 150 bp (~100X coverage).

### 3️⃣ Quality Control  
- **FASTQC**: quality assessment of raw reads.  
- **Prinseq**: trimming of low-quality bases and adapter removal.

### 4️⃣ Genome Assembly  
- **Flye**: assembly of PacBio long reads.  
- **Pilon**: polishing of the assembly using Illumina reads.  
- **QUAST**: assessment of assembly quality (N50, L50, completeness).

### 5️⃣ Genome Annotation  
- **Prokka**: automated annotation of the bacterial genome.  
- **Ensembl Bacteria**: functional annotation and reference data integration.

### 6️⃣ Variant Calling and Analysis  
- **BWA**: alignment of Illumina reads to the reference genome (DSM 428).  
- **Samtools**: manipulation of BAM files (sorting, indexing).  
- **GATK HaplotypeCaller**: SNP and indel calling.  
- **VCFtools**: variant filtration (QUAL>30, DP>10).  
- **SnpEff**: variant effect annotation.  
- **CNV detection**: coverage-based analysis (CNVnator or similar, if required).

### 7️⃣ Visualization and Reporting  
- **IGV (Integrative Genomics Viewer)**: visualization of variant calls and coverage.

---

## Expected Results

- High-quality de novo genome assembly of the industrial strain.  
- Complete list of SNPs, indels and CNVs compared to the reference genome.  
- Functional annotation of variants in PHA biosynthetic genes (phaA, phaB, phaC).  
- Summary tables and visual reports for downstream applications.

---

## Computational Resources

- Linux environment with standard bioinformatics tools (command-line).  
- Sufficient RAM (≥128 GB) for assembly and polishing steps.  
- Storage (~1–2 TB) for raw and processed data.

---

## Budget (Indicative)

| Step | Cost (€) |
|------|----------|
| Sample prep and DNA extraction | 2,500 |
| PacBio sequencing | 50,000 |
| Illumina sequencing | 30,000 |
| Library preparation (PacBio + Illumina) | 10,000 |
| Assembly and variant analysis | 30,000 |
| Personnel (6 months) | 60,000 |
| Dissemination and consumables | 5,000 |
| **Total** | **187,500** |

---

## Contact

For questions, please contact:  
**[Your Name]**  
*Applied Genomics Project – University of Bologna*  
Email: [your-email@example.com]

---



