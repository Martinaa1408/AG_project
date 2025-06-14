# De novo Genome Sequencing and Comparative Genomic Analysis of a Novel Fungal Isolate for Bioeconomic Applications in Bioplastic Degradation

## Table of Contents

- [Project Overview](#project-overview)
- [Background and Rationale](#background-and-rationale)
- [Objectives](#objectives)
- [Methodological Pipeline](#methodological-pipeline)
- [Detailed Cost Breakdown](#detailed-cost-breakdown)
- [Expected Results and Bioeconomy Impact](#expected-results-and-bioeconomy-impact)
- [Sequencing Technology Overview](#sequencing-technology-overview)
- [Data Management](#data-management)
- [References](#references)
- [Contact](#contact)

---

## Project Overview
This project applies advanced genomic techniques to support the circular bioeconomy by identifying fungal genes involved in the enzymatic degradation of biodegradable plastics. A novel fungal isolate was obtained from a compost bin containing PLA (polylactic acid) fragments. Through whole-genome sequencing, de novo assembly, functional annotation, and comparative genomics, this study aims to uncover new enzymes potentially useful in sustainable waste management.

---

## Background and Rationale
The accumulation of plastic waste, even biodegradable plastics like PLA and PHA, represents a growing environmental issue. Fungi are underexplored but promising candidates for biodegradation, thanks to their ability to secrete powerful extracellular enzymes such as esterases and cutinases. This project is aligned with the core themes of Applied Genomics, focusing on sequencing-based genome analysis and functional annotation. By characterizing the genome of a novel fungal strain isolated from a PLA-contaminated site, we aim to contribute to innovative solutions in environmental biotechnology and bioeconomy.

---

## Objectives
- Isolate and characterize a fungal strain from a compostable PLA-contaminated environment.
- Perform high-quality whole-genome sequencing and de novo assembly using hybrid technologies.
- Annotate coding sequences and identify gene families involved in bioplastic degradation.
- Compare the genome to those of known industrial fungal species to highlight unique enzymatic features.
- Provide a publicly available genome resource to support future functional and biotechnological studies.

---

## Methodological Pipeline

### 1. Sample Collection and Fungal Isolation
- Soil and compost samples were collected from a household compost bin containing partially degraded PLA cups.
- Fungal strains growing directly on PLA fragments were isolated on PDA media and screened for esterase-like activity.

### 2. DNA Extraction and QC
- High molecular weight genomic DNA was extracted using a CTAB-based protocol adapted for filamentous fungi.
- DNA quality and concentration were evaluated via Nanodrop, Qubit, and agarose gel electrophoresis.

### 3. Genome Sequencing
- **Short-read sequencing**: Illumina NovaSeq PE150 (~100X coverage) for accurate base-calling.
- **Long-read sequencing**: Oxford Nanopore (MinION/GridION, ~30–50X coverage) to support contiguity and structural resolution.

### 4. Assembly and Quality Assessment
- Long reads assembled using Flye; short reads used to polish the assembly with Pilon.
- Assembly quality assessed via QUAST (N50, L50) and BUSCO with the "fungi_odb10" dataset.

### 5. Genome Annotation
- Gene prediction performed with BRAKER2 and Augustus.
- Functional annotation using InterProScan, eggNOG-mapper, Pfam, and KEGG.
- Identification of putative bioplastic-degrading enzymes based on the presence of relevant domains (e.g., esterase, cutinase, lipase, PHA depolymerase).

### 6. Comparative Genomic Analysis
- Whole-genome comparison with *Trichoderma reesei*, *Aspergillus niger*, and selected Fusarium species.
- Ortholog identification with OrthoFinder and gene family clustering.
- Synteny analysis with MCScanX to highlight conserved and strain-specific regions.

---

## Detailed Cost Breakdown

| Activity | Estimated Cost (€) | Description |
|-------------------------------|--------------------|-------------|
| Sampling & Isolation | 3,000 | Compost handling, fungal plating |
| DNA Extraction & QC | 2,000 | Fungal genomic DNA protocol, reagents |
| Illumina PE150 Sequencing | 30,000 | Short-read sequencing and library prep |
| Oxford Nanopore Sequencing | 50,000 | Long-read sequencing |
| Genome Assembly & Polishing | 10,000 | Hybrid assembly with Flye and Pilon |
| Functional Annotation | 15,000 | BRAKER2, InterProScan, KEGG, etc. |
| Comparative Genomics | 20,000 | Ortholog detection, synteny, tree building |
| Personnel (lab & bioinfo) | 60,000 | Wet and dry lab salaries |
| Reports, Dissemination | 10,000 | Final report, poster, publications |
| **Total** | **200,000** |  |

---

## Expected Results and Bioeconomy Impact
- High-quality annotated genome of a novel fungal isolate from a PLA-contaminated environment.
- Identification of novel enzyme-coding genes involved in plastic degradation.
- Comparative insights into the evolution of fungal biodegradative capacity.
- New genomic data to support the development of fungal-based biodegradation strategies.
- Contribution to European Green Deal objectives by promoting nature-based, genomics-informed waste solutions.

---

## Sequencing Technology Overview

| Technology       | Read Length      | Coverage      | Phred Score | Advantages                                              | Typical Use            |
|------------------|------------------|---------------|-------------|---------------------------------------------------------|------------------------|
| Illumina NovaSeq | ~150 bp (paired) | ~100X         | >30         | High accuracy, ideal for polishing and SNP calling      | Polishing assemblies   |
| Oxford Nanopore  | ~10–20 kb        | ~30–50X       | >30         | Long reads capture repetitive and structural features   | De novo assembly       |

### Additional Concepts (aligned with AG course)
- **Coverage**: ensures redundancy and confidence in variant calling and assembly.
- **Hybrid assembly**: combines the structural accuracy of long reads with the base-level precision of short reads.
- **BUSCO**: benchmark for genome completeness using orthologous genes.
- **Orthology inference**: essential for identifying conserved biodegradation pathways in comparative genomics.
- **Domain-based annotation**: key to enzyme discovery in functional genomics.

---

## Data Management
- Raw sequencing data deposited to ENA under open access BioProject ID.
- Assembly files, annotation tables, and comparative outputs hosted in this GitHub repository.
- All data and scripts managed under FAIR principles and institutional guidelines.

---

## References

- Urbanek, A.K. et al. (2018). Biodegradation of plastics by fungal communities – opportunities and limitations. *Appl Microbiol Biotechnol*.
- Harms, H. et al. (2021). Plastics in the environment – fungal enzymes to the rescue? *Biotechnol Adv*.
- BUSCO: https://busco.ezlab.org/
- BRAKER2: https://github.com/Gaius-Augustus/BRAKER
- InterProScan: https://www.ebi.ac.uk/interpro/
- Applied Genomics materials – University of Bologna

---

## Contact
Martina Castellucci  
Email: martina.castellucci@studio.unibo.it  
Applied Genomics Project – University of Bologna  
