# De novo Genome Sequencing and Comparative Genomic Analysis for Bioeconomic Optimization of Bioplastic Degradation in Fusarium solani

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
This project aims to support sustainable solutions to plastic pollution by exploring the biodegradation potential of fungi. Through the de novo sequencing and comparative analysis of the fungal genome of *Fusarium solani*, isolated from a bioplastic-contaminated soil, the study will identify candidate genes encoding enzymes involved in the degradation of biodegradable plastics such as PLA, PHA, and PCL. This genomic approach contributes to the circular bioeconomy by enabling nature-based waste management strategies.

---

## Background and Rationale
Biodegradable plastics like PLA are increasingly used as eco-friendly alternatives to petroleum-based polymers. However, their natural degradation remains limited without specific microbial or enzymatic interventions. Fungi, particularly *Fusarium solani*, have demonstrated the ability to degrade PLA via extracellular enzymes (e.g., esterases, cutinases). Yet, the genomic determinants of this phenotype remain poorly characterized. By sequencing and analyzing the genome of *F. solani*, this project provides insights into its biotechnological potential for bioplastic degradation.

---

## Objectives
- Perform high-quality de novo sequencing and genome assembly of an environmental *Fusarium solani* strain.
- Annotate the genome with a focus on genes encoding enzymes such as esterases, cutinases, and depolymerases.
- Compare the genome to related fungal species with known enzymatic capabilities (*Trichoderma reesei*, *Aspergillus niger*).
- Identify conserved and unique gene clusters that may contribute to bioplastic degradation.

---

## Methodological Pipeline

### 1. Fungal Isolation and DNA Extraction
- Soil sample collection from a PLA-contaminated environment.
- Fungal isolation via plating and morphological selection.
- Genomic DNA extraction using a fungal-specific protocol (e.g., CTAB-based or commercial kits).
- DNA quality control with Nanodrop, Qubit, and gel electrophoresis.

### 2. Sequencing
- **Illumina NovaSeq (PE150, ~100X coverage)** for high-accuracy short reads.
- **Oxford Nanopore (MinION/GridION, ~30–50X coverage)** for long reads supporting assembly continuity.

### 3. Genome Assembly
- Long-read assembly using Flye.
- Polishing with short reads using Pilon.
- Evaluation of assembly metrics (N50, completeness via BUSCO, total genome size).

### 4. Genome Annotation
- Gene prediction using BRAKER2 and Augustus.
- Functional annotation with InterProScan, eggNOG-mapper, Pfam, and KEGG.
- Targeted search for genes with domains associated with plastic degradation (esterases, cutinases, lipases).

### 5. Comparative Genomics
- Whole-genome comparison against *Trichoderma reesei* and *Aspergillus niger* reference genomes.
- Orthologous cluster identification via OrthoFinder.
- Synteny and gene cluster visualization using MCScanX or similar tools.

---

## Detailed Cost Breakdown

| Activity | Estimated Cost (€) | Description |
|-------------------------------|--------------------|-------------|
| Isolation & DNA Extraction | 5,000 | Soil processing, fungal growth, DNA kits |
| Illumina Sequencing | 30,000 | Paired-end sequencing (NovaSeq) |
| Nanopore Sequencing | 50,000 | Long-read sequencing (MinION/GridION) |
| Assembly and Annotation | 25,000 | Assembly pipelines, annotation, software licenses |
| Comparative Genomics | 20,000 | Reference genome analysis, visualization |
| Personnel | 60,000 | Salary for wet lab and computational staff |
| Dissemination & Consumables | 10,000 | Lab supplies, reporting, data sharing |
| **Total** | **200,000** |  |

---

## Expected Results and Bioeconomy Impact
- A high-quality annotated genome of *Fusarium solani* from a bioplastic environment.
- Identification of candidate genes and pathways for plastic degradation.
- Comparative insights into enzymatic repertoires across relevant fungal species.
- Foundation for applied biotechnology aimed at plastic biodegradation and bioeconomy innovation.

---

## Sequencing Technology Overview

| Technology       | Read Length      | Coverage      | Phred Score | Advantages                                              | Typical Use            |
|------------------|------------------|---------------|-------------|---------------------------------------------------------|------------------------|
| Illumina NovaSeq | ~150 bp (paired) | ~100X         | >30         | Accurate and affordable base calling                    | Polishing and variant calling |
| Oxford Nanopore  | ~15–20 kb        | ~30–50X       | >30         | Captures long genomic repeats; portable                  | De novo assembly       |

---

## Data Management
- Raw reads and final assemblies deposited in ENA/NCBI SRA under open-access BioProject ID.
- Annotated genome, metadata, and analysis scripts maintained in this GitHub repository.
- All data handled under FAIR principles and institutional data protection policies.

---

## References

- Urbanek, A.K. et al. (2018). Biodegradation of plastics by fungal communities – opportunities and limitations. *Applied Microbiology and Biotechnology*.
- Ma, Y. et al. (2021). Biodegradation of polylactic acid by *Fusarium* species: gene and enzyme characterization. *Frontiers in Microbiology*.
- Nanopore sequencing: https://nanoporetech.com  
- InterProScan: https://www.ebi.ac.uk/interpro/
- Applied Genomics course materials – University of Bologna

---

## Contact
Martina Castellucci  
Email: martina.castellucci@studio.unibo.it  
Applied Genomics Project – University of Bologna
