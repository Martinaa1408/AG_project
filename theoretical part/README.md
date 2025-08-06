## SLIDE

- [🔎Google Drive](https://drive.google.com/file/d/1tV58Ldxbase2jpN9sgMsIyVQvUfpziG5/view?usp=sharing)

### 1. Genetics Foundations

Branches of Genetics:

Classical / Transmission genetics → Mendel, crosses, law of segregation & independent assortment.

Molecular genetics → Gene structure, replication, transcription, translation, regulation.

Population genetics → Allele frequencies, Hardy-Weinberg equilibrium, evolution.

Quantitative genetics → Polygenic traits, heritability, variance decomposition.

Key concepts:

Gene, allele, genotype, phenotype.

Pedigree symbols & PLINK .fam file format.

Chromosomes, karyotypes, sex determination systems.

Morgan’s linkage and crossing-over → recombination & mapping.

---

### 2. Genomics & NGS Basics

Genome definitions & history:

Genome = all DNA of an organism.

Genomics = analysis of full gene sets & their function.

Origin: HGP → comparative & functional genomics.

NGS Generations:

Sanger sequencing (dideoxy, 1st-gen).

NGS / 2nd-gen: Illumina, 454, SOLiD, Ion Torrent.

3rd-gen / TGS: PacBio SMRT, Oxford Nanopore.

NGS Key files & metrics:

FASTQ (raw reads), BAM (alignments), VCF (variants), BED (features).

Phred score (Q): Q30 = 0.1% error.

Depth (coverage) = LN/G; Breadth of coverage %.

---

### 3. NGS Data Analysis Pipeline

Typical workflow:

QC → FastQC, trimming (PrinSeq, Trimmomatic).

Alignment → BWA, Bowtie; BAM/SAM; mapping quality (MQ).

Filtering → Remove duplicates (Picard), low MQ reads.

Variant calling → GATK, VarScan, bcftools; output VCF.

Variant annotation → SnpEff, ANNOVAR, Ensembl VEP.

Manual inspection → IGV.

QC checks:

Base quality per cycle.

GC content distribution.

Read duplication levels.

---

### 4. Genome Assembly

Strategies:

De novo vs Reference-guided.

Shotgun approaches: Classical, hierarchical BAC-based.

Graph models:

OLC (Overlap-Layout-Consensus)

De Bruijn Graphs (k-mer nodes/edges, Eulerian path).

Hybrid assemblies (Illumina + TGS).

Quality metrics:

N50, NG50, L50 (contiguity).

BUSCO → genome completeness.

Gap filling & scaffolding → paired-end, mate-pair, Hi-C, optical maps.

---

### 5. Genome Annotation

Steps:

Repeat masking → RepeatMasker, DFAM, REPET.

Structural annotation → ab-initio (AUGUSTUS), extrinsic (RNA-seq, protein homology), combiners.

Functional annotation → GO terms, BLAST, UniProt, domains.

Output formats:

GFF3, BED, GenBank/EMBL for genome browsers.

---

### 6. Transcriptomics & RNA-Seq

Principles:

mRNA enrichment → poly-A selection.

RNA → cDNA → library prep → NGS.

Reads can be SE or PE; align-then-assemble vs de novo assembly.

Expression analysis:

Normalization metrics: RPKM / FPKM / TPM.

Differential expression: DESeq2, edgeR.

Junction-aware mapping → STAR, HISAT2.

---

### 7. Population Genomics & GWAS

PLINK file types:

.ped + .map (text); .bed/.bim/.fam (binary).

Quality control filters: --mind, --geno, --maf, --hwe.

Population analysis:

MAF, HWE, FST.

PCA / MDS for population structure.

GWAS & ROH (Runs of Homozygosity).

Visualization: Manhattan plots, QQ plots.

---

### 8. High-Throughput Genotyping & CNV

Technologies:

SNP arrays (Illumina BeadChip, Axiom arrays).

NGS-based genotyping: RAD-Seq, ddRAD, GBS.

Copy Number Variation (CNV) detection: aCGH, SNP intensity, NGS depth.

---

### 9. Epigenomics & Functional Genomics

ChIP-Seq → protein-DNA binding sites.

Methyl-Seq / WGBS → CpG methylation & regulation.

ATAC-seq / DNase-seq → open chromatin regions.

AntiSMASH / functional pipelines → BGCs and secondary metabolites.

---

### 10. Practical NGS Case Study

Pipeline with Galaxy:

FASTQ → QC → Trimming → BWA → BAM → Dedup → Variant Calling → VCF → Annotation with VEP.

Key Takeaways:

Always collect sample & library metadata.

Evaluate sequencing depth & breadth.

Validate candidate SNPs/indels with visualization (IGV).

