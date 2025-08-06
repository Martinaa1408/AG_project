# Applied Genomics – Slide Recap

Full slide deck: [Google Drive](https://drive.google.com/file/d/1tV58Ldxbase2jpN9sgMsIyVQvUfpziG5/view?usp=sharing)

This document summarizes the **core concepts of the Applied Genomics course**, including classical and molecular genetics, NGS technologies, genome assembly, functional and population genomics, and applied bioinformatics workflows.

---

## 1️⃣ Genetics Foundations

- **Branches of Genetics**
  - **Classical / Transmission** → Mendel’s laws, segregation, independent assortment
  - **Molecular** → Gene structure, replication, transcription, translation, regulation
  - **Population** → Allele frequencies, Hardy–Weinberg equilibrium (HWE), evolution
  - **Quantitative** → Polygenic traits, heritability, variance decomposition

- **Key Concepts**
  - Gene / allele / genotype / phenotype  
  - Pedigree symbols → male ◻, female ●, carriers half-filled  
  - PLINK `.fam` format → Family ID, individual ID, parental info, sex, phenotype
  - Morgan’s linkage → crossing-over & genetic mapping

---

## 2️⃣ Genomics & NGS Essentials

- **Genome & Genomics**
  - Genome = complete DNA content  
  - Genomics = analysis of all genes and functions  
  - Evolution from **HGP** → comparative & functional genomics

- **Sequencing Generations**
  - 1st: **Sanger** (dideoxy chain-termination)  
  - 2nd: **Illumina, 454, SOLiD, Ion Torrent**  
  - 3rd: **PacBio SMRT, Oxford Nanopore (ONT)**

- **Core NGS Files & Metrics**
  - **FASTQ** (raw reads + Phred quality)  
  - **BAM/SAM** (aligned reads)  
  - **VCF** (variants), **BED** (features)  
  - Phred Q30 ≈ 0.1% error; Coverage = (read length × reads) / genome size

---

## 3️⃣ NGS Data Analysis Pipeline

1. **Quality control** → FastQC, MultiQC  
2. **Trimming / filtering** → Trimmomatic, PrinSeq  
3. **Alignment** → BWA / Bowtie2 → BAM/SAM  
4. **Post-processing** → Sort, index, mark duplicates  
5. **Variant calling** → GATK, FreeBayes, VarScan → VCF  
6. **Annotation** → SnpEff, ANNOVAR, Ensembl VEP  
7. **Manual review** → IGV genome browser

**QC checks:** base quality per cycle, GC bias, duplication levels

---

## 4️⃣ Genome Assembly

- **Approaches**
  - **De novo** → no reference genome  
  - **Reference-guided** → use existing genome  
  - **Hybrid** → combine Illumina + ONT/PacBio

- **Graph Models**
  - OLC (Overlap–Layout–Consensus)  
  - **De Bruijn Graph** → (k-1)-mers as nodes, k-mers as edges

- **Key Metrics**
  - **N50 / NG50 / L50** → contiguity  
  - **BUSCO** → completeness using single-copy orthologs  
  - **Scaffolding** → paired-end, mate-pair, Hi-C, optical maps

---

## 5️⃣ Genome Annotation

- **Repeat masking** → RepeatMasker, DFAM  
- **Structural annotation** → AUGUSTUS, GeneMark, RNA-Seq evidence  
- **Functional annotation** → BLAST, UniProt, InterProScan, GO terms  
- **Outputs** → GFF3, BED, GenBank for browsers

---

## 6️⃣ Transcriptomics & RNA-Seq

- **Workflow**
  - Poly-A mRNA enrichment or rRNA depletion  
  - cDNA library prep → SE/PE sequencing  
  - Alignment → STAR, HISAT2 (splice-aware)  
  - Counting → HTSeq, featureCounts  
  - Differential expression → DESeq2, edgeR

- **Normalization metrics** → RPKM / FPKM / TPM

- **Challenges** → splice junctions, batch effects, low expression genes

---

## 7️⃣ Population Genomics & GWAS

- **PLINK basics**
  - Text: `.ped + .map`  
  - Binary: `.bed + .bim + .fam`  
  - QC filters: `--mind` (missing), `--geno`, `--maf`, `--hwe`

- **Population analyses**
  - MAF, HWE, \( F_{ST} \)  
  - PCA / MDS → population structure

- **GWAS & ROH**
  - Runs of Homozygosity → inbreeding detection  
  - Manhattan and QQ plots for association signals

---

## 8️⃣ High-Throughput Genotyping & CNV

- **SNP arrays** → Illumina BeadChip, Axiom arrays  
- **NGS-based genotyping** → RAD-Seq, ddRAD, GBS  
- **CNV detection**
  - Array CGH (aCGH) → hybridization intensity  
  - Read depth-based → CNVnator, XHMM

---

## 9️⃣ Epigenomics & Functional Genomics

- **ChIP-Seq** → Protein-DNA interactions (TFs, histones)  
- **Methyl-Seq / WGBS** → CpG methylation & epigenetic regulation  
- **ATAC-Seq / DNase-Seq** → Chromatin accessibility  
- **AntiSMASH** → Biosynthetic gene cluster (BGC) prediction

---

## 🔟 Practical NGS Case Study

**Pipeline Example (Galaxy / Snakemake)**  

FASTQ → QC → Trimming → Alignment (BWA)
→ BAM processing → Variant Calling (GATK)
→ VCF → Annotation (VEP) → Visualization (IGV)


**Key Takeaways**
- Always collect **metadata** (sample, library, environment)  
- Evaluate **depth & breadth** before interpretation  
- Validate **candidate variants** visually in IGV  

