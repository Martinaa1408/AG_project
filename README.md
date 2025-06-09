## De novo Genome Sequencing and Comparative Genomic Analysis of *Cupriavidus necator* H16 for Bioeconomy Applications


## Table of Contents

- [Project Overview](#project-overview)
- [Background and Rationale](#background-and-rationale)
- [Objectives](#objectives)
- [Methodological Pipeline](#methodological-pipeline)
- [Detailed Cost Breakdown](#detailed-cost-breakdown)
- [Expected Results and Bioeconomy Impact](#expected-results-and-bioeconomy-impact)
- [Data Management](#data-management)
- [Contact](#contact)

---

## Project Overview
This project focuses on the **de novo genome sequencing** of an industrial strain of *Cupriavidus necator* H16, combined with a **comparative analysis** against the reference strain (DSM 428) and selected industrial strains. The aim is to identify key genomic features relevant to sustainable bioplastic production, aligning with the principles of the bioeconomy.

---

## Background and Rationale
*Cupriavidus necator* H16 is an industrially relevant bacterium known for producing polyhydroxyalkanoates (PHAs), a family of biodegradable plastics. The genome of the reference strain (DSM 428) is publicly available, but industrial strains may present genomic variations that influence PHA yields and stress tolerance. Comparative genomic analysis is crucial to identify key differences, guiding future metabolic engineering and contributing to circular economy initiatives.

This approach ensures that the project is not solely focused on bacteria but integrates **genome assembly, annotation, and comparative genomics**, thus aligning with the recommended project guidelines.

---

## Objectives
- Assemble a high-quality genome of an industrial strain of *C. necator* H16 using long-read and short-read sequencing.
- Annotate the genome, focusing on genes involved in PHA production and relevant metabolic pathways.
- Compare the assembled genome with the reference strain DSM 428 and with selected industrial strains to identify genomic differences that may impact bioplastic production.
- Contribute to sustainable industrial processes through applied genomics research.

---

## Methodological Pipeline

### 1. Sample Collection and DNA Extraction
- Culture the industrial strain in appropriate media under controlled conditions.
- Extract high-molecular-weight DNA using commercial kits (e.g., Qiagen DNeasy).
- Assess DNA purity (Nanodrop) and concentration (Qubit).

### 2. Sequencing
- **Short-read sequencing**: Illumina NovaSeq PE150 (~100X coverage) for high-accuracy base calls.
- **Long-read sequencing**: PacBio HiFi (~50-100X coverage) to improve assembly contiguity.

### 3. Genome Assembly
- Assemble long reads with Flye or a similar assembler.
- Polish the assembly using Illumina reads with Pilon.
- Evaluate assembly quality using QUAST (N50, L50, completeness).

### 4. Genome Annotation
- Annotate genes using Prokka.
- Cross-reference PHA-related genes (e.g., phaA, phaB, phaC) with public databases and literature.

### 5. Comparative Genomics
- Compare the assembled genome with:
  - The DSM 428 reference genome.
  - Additional industrial strains available in public repositories (e.g., NCBI).
- Use tools like Mauve or MUMmer to identify synteny, insertions/deletions, and gene cluster differences.

---

## Detailed Cost Breakdown

| Activity | Estimated Cost (€) | Description |
|-------------------------------|--------------------|-------------|
| Sample Preparation & DNA Extraction | 5,000 | Media, DNA extraction kits, quality checks (Nanodrop, Qubit). |
| Illumina Sequencing | 35,000 | Short-read sequencing with library preparation. |
| PacBio Sequencing | 50,000 | Long-read sequencing with library preparation. |
| Assembly and Annotation | 20,000 | Genome assembly, polishing, and annotation with essential software. |
| Comparative Analysis | 15,000 | Identification of gene differences, synteny, and visualization. |
| Personnel (essential) | 60,000 | Researcher(s) dedicated to lab work, data curation, and report preparation. |
| Dissemination & Consumables | 15,000 | Reports, presentations, laboratory supplies. |
| **Total** | **200,000** |  |

---

## Expected Results and Bioeconomy Impact
- Generation of a complete, high-quality genome assembly of the industrial strain.
- Identification of genes and pathways relevant to PHA production and stress adaptation.
- Comparative insights that support strain optimization and sustainable bioplastic production.
- Contribution to the European bioeconomy strategy through applied genomics research.

---

## Data Management
- All sequencing data and assemblies will be stored securely on institutional servers.
- Data will be submitted to NCBI SRA/ENA for public access.
- Reports, figures, and scripts will be stored in a dedicated GitHub repository.

---

## Contact
Martina Castellucci 

Applied Genomics Project – University of Bologna  
Email: martina.castellucci@studio.unibo.it
