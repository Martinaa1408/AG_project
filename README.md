# De novo Genome Sequencing, Transcriptomics, and Functional Annotation of a PLA-Degrading Fungal Isolate for Bioeconomic Innovation

## Table of Contents

* [Project Overview](#project-overview)
* [Background and Rationale](#background-and-rationale)
* [Objectives](#objectives)
* [Methodological Pipeline](#methodological-pipeline)
* [Detailed Cost Breakdown](#detailed-cost-breakdown)
* [Expected Results and Bioeconomy Impact](#expected-results-and-bioeconomy-impact)
* [Sequencing Technology Overview](#sequencing-technology-overview)
* [Data Management](#data-management)
* [References](#references)
* [Contact](#contact)

---

## Project Overview

This project applies advanced genomic and transcriptomic techniques to support the circular bioeconomy by identifying fungal genes involved in the enzymatic degradation of biodegradable plastics. A novel fungal isolate was obtained from a compost bin containing PLA (polylactic acid) fragments. Through whole-genome sequencing, de novo assembly, transcriptomic profiling, functional annotation, and comparative genomics, this study aims to uncover new enzymes potentially useful in sustainable waste management.

---

## Background and Rationale

The accumulation of plastic waste, even biodegradable plastics like PLA and PHA, represents a growing environmental issue. Fungi are underexplored but promising candidates for biodegradation, thanks to their ability to secrete powerful extracellular enzymes such as esterases and cutinases. This project is aligned with the core themes of Applied Genomics, focusing on sequencing-based genome analysis, transcriptomics, and functional annotation. By characterizing the genome and gene expression of a novel fungal strain isolated from a PLA-contaminated site, we aim to contribute to innovative solutions in environmental biotechnology and bioeconomy.

---

## Objectives

* Isolate and characterize a fungal strain from a compostable PLA-contaminated environment.
* Perform high-quality whole-genome sequencing and de novo assembly using hybrid technologies.
* Annotate coding sequences and identify gene families involved in bioplastic degradation.
* Perform transcriptomic profiling under PLA/PHA-containing media vs control to detect differentially expressed genes.
* Compare the genome to those of known industrial fungal species and non-biodegrading relatives to highlight lineage-specific features.
* Identify biosynthetic gene clusters (BGCs) and lineage-specific metabolic pathways using AntiSMASH and MAUVE.
* Provide a publicly available genome and transcriptome resource to support future biotechnological studies.

---

## Methodological Pipeline

### 1. Sample Collection and Fungal Isolation

* Soil and compost samples were collected from a household compost bin containing partially degraded PLA cups.
* Fungal strains growing directly on PLA fragments were isolated on PDA media and screened for esterase-like activity.

### 2. DNA Extraction and QC

* High molecular weight genomic DNA was extracted using a CTAB-based protocol adapted for filamentous fungi.
* DNA quality and concentration were evaluated via Nanodrop, Qubit, and agarose gel electrophoresis.

### 3. Genome Sequencing

* **Short-read sequencing**: Illumina NovaSeq PE150 (\~100X coverage) for accurate base-calling.
* **Long-read sequencing**: Oxford Nanopore (MinION/GridION, \~30–50X coverage) to support contiguity and structural resolution.

### 4. Assembly and Quality Assessment

* Long reads assembled using Flye; short reads used to polish the assembly with Pilon.
* Assembly quality assessed via QUAST (N50, L50) and BUSCO with the "fungi\_odb10" dataset.

### 5. Transcriptome Analysis

* RNA was extracted from fungal cultures grown in standard media and PLA/PHA-supplemented media.
* RNA-Seq libraries prepared and sequenced with Illumina.
* Reads aligned to genome using STAR, and expression levels quantified with featureCounts.
* Differential expression analysis using DESeq2 to identify overexpressed genes in bioplastic conditions.

### 6. Genome Annotation

* Gene prediction performed with **MAKER3**, integrating evidence from transcriptomic reads, Augustus, and GeneMark.
* Functional annotation using InterProScan, eggNOG-mapper, Pfam, KEGG, and **Fungal AntiSMASH** for BGCs.
* Identification of candidate bioplastic-degrading enzymes based on relevant domains (e.g., esterase, cutinase, lipase, PHA depolymerase), refined with **Markov models** for enzyme families.

### 7. Comparative Genomic Analysis

* Whole-genome comparison with *Trichoderma reesei*, *Aspergillus niger*, selected *Fusarium* species, and non-biodegrading close relatives.
* Ortholog identification with OrthoFinder and gene family clustering.
* Synteny analysis and alignment with **MAUVE** (at varying window sizes) and **MCScanX** to highlight lineage-specific genomic regions.

### 8. Phylogenetic Analysis

* Single-copy orthologs identified with OrthoFinder.
* Aligned with MAFFT, conserved blocks selected via Gblocks, concatenated using AMAS.
* Phylogenetic trees constructed using MEGA11, FastTree, RAxML, and MrBayes.

---

## Detailed Cost Breakdown

| Personnel              | Activity                               | Estimated Cost (€) | Description                                       |
| ---------------------- | -------------------------------------- | ------------------ | ------------------------------------------------- |
| Wet lab Postdoc salary |                                        | 45,000             | One year project                                  |
|                        | Sampling & Isolation                   | 3,000              | Compost handling, fungal plating                  |
|                        | DNA Extraction & QC                    | 2,000              | Fungal genomic DNA protocol, reagents             |
|                        | Illumina PE150 Sequencing              | 30,000             | Short-read sequencing and library prep            |
|                        | Oxford Nanopore Sequencing             | 50,000             | Long-read sequencing                              |
|                        | RNA extraction (12 libraries)          | 1,000              |                                                   |
|                        | RNA library preparation (12 libraries) | 1,500              |                                                   |
|                        | RNASeq Illumina (12 libraries)         | 2,500              |                                                   |
| Bioinformatics Postdoc |                                        | 45,000             | One year project                                  |
|                        | Genome Assembly & Polishing            | 0                  | Hybrid assembly with Flye and Pilon               |
|                        | Functional Annotation                  | 0                  | MAKER3, InterProScan, KEGG, AntiSMASH             |
|                        | Comparative Genomics                   | 0                  | Orthologs, synteny, MAUVE, lineage-specific genes |
|                        | Non-open source softwares              | 5,000              | Genious & others                                  |
|                        | Reports, Dissemination                 | 15,000             | Publications & Conferences                        |
|                        | **Total**                              | **200,000**        |                                                   |

---

## Expected Results and Bioeconomy Impact

* High-quality annotated genome and transcriptome of a novel fungal isolate from a PLA-contaminated environment.
* Identification and expression analysis of novel enzyme-coding genes involved in plastic degradation.
* Detection of biosynthetic gene clusters (BGCs) with potential metabolic innovation.
* Comparative insights into lineage-specific genomic features related to biodegradation.
* Phylogenetic positioning of the isolate using robust single-copy ortholog pipelines.
* Contribution to European Green Deal objectives through innovative genomic tools for waste management.

---

## Sequencing Technology Overview

| Technology       | Read Length       | Coverage | Phred Score | Advantages                                            | Typical Use          |
| ---------------- | ----------------- | -------- | ----------- | ----------------------------------------------------- | -------------------- |
| Illumina NovaSeq | \~150 bp (paired) | \~100X   | >30         | High accuracy, ideal for polishing and SNP calling    | Polishing assemblies |
| Oxford Nanopore  | \~10–20 kb        | \~30–50X | >30         | Long reads capture repetitive and structural features | De novo assembly     |

### Additional Concepts (aligned with AG course)

* **Coverage**: ensures redundancy and confidence in variant calling and assembly.
* **Hybrid assembly**: combines the structural accuracy of long reads with the base-level precision of short reads.
* **BUSCO**: benchmark for genome completeness using orthologous genes.
* **Orthology inference**: essential for identifying conserved biodegradation pathways in comparative genomics.
* **Domain-based annotation**: key to enzyme discovery in functional genomics.
* **Markov models**: used to refine identification of enzyme families from domain structures.
* **MAKER3**: genome annotation framework integrating transcript and ab initio predictions.

---

## Data Management

* Raw sequencing and RNA-Seq data deposited to ENA under open access BioProject ID.
* Assembly files, annotation tables, expression matrices, and comparative outputs hosted in this GitHub repository.
* All data and scripts managed under FAIR principles and institutional guidelines.

---

## References

* Urbanek, A.K. et al. (2018). Biodegradation of plastics by fungal communities – opportunities and limitations. *Appl Microbiol Biotechnol*.
* Harms, H. et al. (2021). Plastics in the environment – fungal enzymes to the rescue? *Biotechnol Adv*.
* BUSCO: [https://busco.ezlab.org/](https://busco.ezlab.org/)
* BRAKER2: [https://github.com/Gaius-Augustus/BRAKER](https://github.com/Gaius-Augustus/BRAKER)
* InterProScan: [https://www.ebi.ac.uk/interpro/](https://www.ebi.ac.uk/interpro/)
* AntiSMASH: [https://fungismash.secondarymetabolites.org/](https://fungismash.secondarymetabolites.org/)
* MAKER3: [http://www.yandell-lab.org/software/maker.html](http://www.yandell-lab.org/software/maker.html)
* Applied Genomics materials – University of Bologna

---

## Contact

Martina Castellucci
Email: [martina.castellucci@studio.unibo.it](mailto:martina.castellucci@studio.unibo.it)
Applied Genomics Project – University of Bologna
