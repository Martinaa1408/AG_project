## De novo Genome Sequencing and Comparative Genomic Analysis for Bioeconomic Optimization of Bioplastic Production


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
This project supports the sustainable production of biodegradable plastics (polyhydroxyalkanoates, PHAs) through applied genomics. By sequencing and comparing the genome of an industrial strain to the reference strain (DSM 428) and additional industrial strains, the project aims to identify key genes and pathways relevant to industrial PHA production, aligning with circular bioeconomy goals.

### About *Cupriavidus necator* H16

*Cupriavidus necator* H16 (DSM 428) is a Gram-negative, facultative chemolithoautotrophic bacterium known for its remarkable metabolic versatility. It can grow on organic substrates (heterotrophy) or by fixing carbon dioxide using hydrogen as an energy source (autotrophy). Notably, C. necator H16 is widely recognized for its ability to produce polyhydroxyalkanoates (PHAs), biodegradable bioplastics of significant interest for the bioeconomy. Its genome comprises two chromosomes and a megaplasmid, with a total size of approximately 7.4 Mb. The well-characterized genome makes it an excellent model for studies on industrial bioproduction of bioplastics and sustainable bioprocesses.

---

## Background and Rationale
PHAs are industrially relevant bioplastics that reduce reliance on fossil fuels. Industrial production using bacterial fermentation is an established process, yet strain performance can vary significantly. Integrating **de novo genome assembly, annotation, and comparative analysis** allows the identification of genomic differences that may affect PHA yields and stress tolerance, supporting targeted strain improvement strategies aligned with bioeconomic sustainability.

---

## Objectives
- Generate a high-quality genome assembly of an industrial strain used in PHA production using long- and short-read sequencing.
- Annotate the genome with a focus on PHA biosynthesis pathways and stress-response genes.
- Compare the industrial strain’s genome to the DSM 428 reference and other industrial strains to identify relevant genomic differences.
- Apply findings to support sustainable and efficient bioplastic production.

---

## Methodological Pipeline

### 1. Sample Preparation and DNA Extraction
- Grow the industrial strain in standard laboratory conditions.
- Extract high-quality genomic DNA using Qiagen DNeasy or similar kits.
- Assess DNA quality and concentration using Nanodrop and Qubit.

### 2. Sequencing
- **Short-read sequencing**: Illumina NovaSeq PE150 (~100X coverage) for accurate base calling.
- **Long-read sequencing**: PacBio HiFi (~50-100X coverage) to enhance assembly quality.

### 3. Genome Assembly
- Assemble long reads using a suitable assembler (e.g. Flye or equivalent).
- Polish the assembly using Illumina reads with Pilon.
- Assess assembly quality using QUAST (N50, L50, completeness).

### 4. Genome Annotation
- Annotate the genome using Prokka.
- Focus on PHA-related genes (phaA, phaB, phaC) and stress response genes.

### 5. Comparative Genomics
- Compare the assembled genome to:
  - The DSM 428 reference genome.
  - Other industrial strains from public repositories (e.g. NCBI).
- Identify gene differences, insertions/deletions, and synteny using alignment tools.

---

## Detailed Cost Breakdown

| Activity | Estimated Cost (€) | Description |
|-------------------------------|--------------------|-------------|
| Sample Preparation & DNA Extraction | 5,000 | Media, extraction kits, quality checks. |
| Illumina Sequencing | 35,000 | Short-read sequencing and library preparation. |
| PacBio Sequencing | 50,000 | Long-read sequencing and library preparation. |
| Assembly and Annotation | 20,000 | Genome assembly, polishing, and annotation. |
| Comparative Analysis | 15,000 | Gene cluster analysis, synteny, visualization. |
| Personnel (essential) | 60,000 | Lab work, data management, report preparation. |
| Dissemination & Consumables | 15,000 | Reports, presentations, lab supplies. |
| **Total** | **200,000** |  |

---

## Expected Results and Bioeconomy Impact
- High-quality genome assembly of the industrial strain.
- Annotation of PHA biosynthesis genes and stress-related pathways.
- Comparative insights to support strain optimization and sustainable production.
- Contribution to the European bioeconomy strategy by enabling improved bioplastic production.

---

## Sequencing Technology Overview

| Technology       | Read Length      | Coverage      | Phred Score | Advantages                                              | Typical Use            |
|------------------|------------------|---------------|-------------|---------------------------------------------------------|------------------------|
| Illumina NovaSeq | ~150 bp (paired) | ~100X         | >30         | Accurate base calling; cost-effective; high throughput  | Variant calling; polishing |
| PacBio HiFi      | ~15–20 kb        | ~50–100X      | >30         | Long reads resolve repeats; high consensus accuracy     | De novo assembly       |

### Additional Concepts
- **Coverage**: higher coverage increases confidence in genome assembly and variant calling.  
- **Phred Score**: indicates sequencing quality; >30 means 99.9% accuracy.  
- **Hybrid Assembly**: combines long reads for structure and short reads for polishing.  
- **Annotation**: identifies coding sequences, regulatory elements, and pathways (e.g., phaA, phaB, phaC).  
- **Comparative Genomics**: highlights differences and similarities with DSM 428 and industrial strains, aiding in strain optimization.


---

## Data Management
- Data stored on secure institutional servers.
- Raw reads and assemblies deposited in public repositories (e.g. NCBI SRA/ENA).
- Reports, figures, and protocols managed in a dedicated GitHub repository.

---

## References

- Lu, Y. et al. (2019). Complete genome sequence of *Cupriavidus necator* H16 (DSM 428). Microbiology Resource Announcements.  
  https://journals.asm.org/doi/10.1128/mra.00814-19
- *Cupriavidus necator* H16 (DSM 428) — DSMZ Strain Information:  
  https://bacmedia.dsmz.de/strains/view/DSM%20428
- Applied Genomics course materials – University of Bologna.
- Illumina NovaSeq: https://www.illumina.com/systems/sequencing-platforms/novaseq.html
- PacBio HiFi sequencing: https://www.pacb.com/hifi-sequencing/
- Prokka: https://github.com/tseemann/prokka
- BWA: http://bio-bwa.sourceforge.net/
- Samtools: http://www.htslib.org/
- GATK HaplotypeCaller: https://gatk.broadinstitute.org/
- VCFtools: https://vcftools.github.io/
- FastQC: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/

---

## Contact
Martina Castellucci  
Email: martina.castellucci@studio.unibo.it
Applied Genomics Project – University of Bologna  
