## De novo Genome Assembly and Variant Analysis of *Cupriavidus necator* H16


## Table of Contents

- [Overview](#overview)
- [Objectives](#objectives)
- [Materials and Methods](#materials-and-methods)
- [Expected Results](#expected-results)
- [Budget and Explanations](#budget-and-explanations)
- [Additional Notes on Budget Items](#additional-notes-on-budget-items)
- [References](#references)
- [Contact](#contact)

---

## Overview
This project aims to generate a high-quality de novo genome assembly of an industrial strain of *Cupriavidus necator* H16, a key model organism for the production of biodegradable plastics (polyhydroxyalkanoates, PHA). The genome of the DSM 428 reference strain is already available (Lu et al., 2019), but industrial strains may exhibit genetic differences impacting PHA productivity. By sequencing, assembling, and analyzing the genome, we aim to identify mutations (SNPs and CNVs) in genes like phaA, phaB, and phaC that could be exploited for strain improvement.

---

## Objectives
- Assemble a high-quality genome of an industrial strain of *C. necator* H16 using long-read and short-read sequencing.
- Compare the assembled genome with the DSM 428 reference genome to identify SNPs, indels, and CNVs.
- Annotate variants in key PHA biosynthesis genes to assess their potential impact on PHA production.

---

## Materials and Methods

### Sample Collection and DNA Extraction
- An industrial strain of *C. necator* H16 will be grown in LB medium at 30°C.
- Genomic DNA will be extracted using Qiagen DNeasy Blood & Tissue Kit.
- Quality of the DNA will be assessed using Nanodrop (A260/280), Qubit, and agarose gel electrophoresis.

### Sequencing
- Long-read sequencing: PacBio HiFi (~100X genome coverage).
- Short-read sequencing: Illumina NovaSeq 6000 PE150 (~100X genome coverage).

### Quality Control
- FastQC: assessment of raw read quality.
- Prinseq: trimming of low-quality reads and adapters.

### Genome Assembly
- Flye: assembly of long reads.
- Pilon: polishing of the assembly with Illumina short reads.
- QUAST: assessment of assembly quality (N50, L50, completeness).

### Genome Annotation
- Prokka: rapid bacterial genome annotation.
- Ensembl Bacteria: functional reference and confirmation of key genes.

### Variant Analysis
- BWA: alignment of Illumina reads to the reference genome DSM 428.
- Samtools: manipulation of BAM files (sorting, indexing).
- GATK HaplotypeCaller: variant calling for SNPs and indels.
- VCFtools: filtering of variant calls (QUAL >30, DP >10).
- SnpEff: annotation of variant effects.
- CNVnator (or coverage-based analysis): detection of CNVs.

### Visualization
- IGV: visualization of variant calls and read alignments.

---

## Expected Results
- A high-quality genome assembly of the industrial strain of *C. necator* H16.
- A list of SNPs, indels, and CNVs compared to the DSM 428 reference.
- Annotated variants in PHA biosynthesis genes (phaA, phaB, phaC).
- Tables and summary figures to support downstream analyses.

---

## Budget and Explanations

| Step                               | Estimated Cost (€) | Description |
|------------------------------------|--------------------|-------------|
| Sample prep & DNA extraction       | 2,000              | Includes culture media, DNA extraction kits, quality control reagents (Nanodrop, Qubit) and laboratory consumables. |
| PacBio sequencing                  | 60,000             | Covers library preparation (SMRTbell kit) and sequencing (~100X genome coverage). Essential for long-read assembly. |
| Illumina sequencing                | 30,000             | Includes library preparation and sequencing (PE150, ~100X coverage). Used for polishing and accurate variant detection. |
| Library preparation kits           | 10,000             | Specific kits for both PacBio and Illumina library construction. |
| Assembly and variant analysis      | 30,000             | Covers bioinformatics software licensing, computational time, storage, and personnel time for genome assembly, variant calling, and annotation. |
| Personnel (6 months)               | 60,000             | Researcher(s) dedicated to lab work, data analysis, and reporting. |
| Dissemination & consumables        | 8,000              | Costs for report preparation, figures, presentations, and general laboratory consumables (pipettes, gloves, safety equipment). |
| **Total**                          | **200,000**        |  |

---

## Additional Notes on Budget Items
- **Sample prep & DNA extraction**: High-quality DNA is critical for PacBio sequencing success.
- **PacBio sequencing**: Long reads help assemble repetitive regions and reduce assembly fragmentation.
- **Illumina sequencing**: High-accuracy short reads are essential for polishing the assembly and for robust variant calling.
- **Library preparation kits**: Essential for constructing sequencing-ready DNA libraries with appropriate quality and insert size.
- **Assembly and variant analysis**: Includes costs for software tools (some may require licenses or cloud resources) and for trained bioinformaticians.
- **Personnel**: Reflects a half-year contract for one researcher, including social costs and benefits.
- **Dissemination & consumables**: Covers publication fees, figure preparation, basic lab materials, and safety compliance.

---

## References

- Lu, Y. et al. (2019). Complete genome sequence of *Cupriavidus necator* H16 (DSM 428). Microbiology Resource Announcements. https://journals.asm.org/doi/10.1128/mra.00814-19  
- Applied Genomics course notes, University of Bologna.  
- Official documentation of key bioinformatics tools:  
  - Flye: https://github.com/fenderglass/Flye  
  - Pilon: https://github.com/broadinstitute/pilon  
  - BWA: http://bio-bwa.sourceforge.net/  
  - Samtools: http://www.htslib.org/  
  - GATK HaplotypeCaller: https://gatk.broadinstitute.org/  
  - VCFtools: https://vcftools.github.io/  
  - SnpEff: http://snpeff.sourceforge.net/  
  - IGV: https://software.broadinstitute.org/software/igv/  

---

## Contact

Martina Castellucci 

Applied Genomics Project – University of Bologna  
Email: martina.castellucci@studio.unibo.it
