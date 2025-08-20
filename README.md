# Hybrid Genome Assembly and Integrative Multi-Omics of Compost-Derived Fungal Isolates: Insights into Bioplastic Degradation Potential

---

## Table of Contents
- [Scientific Context](#scientific-context)
- [Project Overview](#project-overview)
- [Background and Rationale](#background-and-rationale)
- [Objectives](#objectives)
- [Methodological Pipeline](#methodological-pipeline)
- [Summary of Experimental Workflow](#summary-of-experimental-workflow)
- [Detailed Cost Breakdown](#detailed-cost-breakdown)
- [Expected Key Results](#expected-key-results)
- [Fungal Candidates in PLA/PHA Waste (Simulated Preview)](#fungal-candidates-in-plapha-waste-simulated-preview)
- [Bioeconomy Impact and EEA Policy Alignment](#bioeconomy-impact-and-eea-policy-alignment)
- [Sequencing Technology Overview](#sequencing-technology-overview)
- [References](#references)
- [Contact](#contact)

---

## Scientific Context

Bioplastics such as **PLA (polylactic acid)** and **PHA (polyhydroxyalkanoates)** are sustainable alternatives to petroleum plastics but degrade slowly under natural composting conditions.  
Filamentous fungi are strong candidates for bioplastic degradation because they secrete **hydrolytic enzymes** (esterases, cutinases, lipases) capable of breaking ester bonds.  

PLA-enriched compost provides a **hotspot** for isolating potential degradative fungi.  
This project explores fungal diversity from compost, aiming to identify isolates and genes linked to **bioplastic depolymerization**.

---

## Project Overview

- Compost sampling yielded ~20 fungal isolates.  
- Screening with PDA + PLA + Rhodamine B identified ~5 isolates with clear degradative activity.  
- Selected isolates undergo **multi-omics analysis** (genome, transcriptome, comparative genomics) to reveal candidate enzymes for plastic degradation.  

---

## Background and Rationale

- Bioplastics are increasingly used but not always fully degraded under composting conditions.  
- Fungi adapt well to compost niches and secrete extracellular enzymes capable of attacking polymers.  
- A lack of **high-quality fungal genomes** and **integrated transcriptomic data** hinders enzyme discovery.  
- This project uses **hybrid sequencing and integrative multi-omics** to bridge this gap.  

---

## Objectives

- Isolate and identify fungal strains from PLA/PHA compost.  
- Generate **hybrid genome assemblies** (Illumina + Nanopore).  
- Perform **functional annotation** focusing on CAZymes and degradative enzymes.  
- Conduct **comparative genomics** between degraders and non-degraders.  
- Apply **RNA-Seq** to detect genes induced under PLA/PHA.  
- Integrate datasets to provide a **pipeline for enzyme discovery**.  

---

## Methodological Pipeline

1. **Sample Collection & Isolation** – Compost plating on PDA + PLA + Rhodamine B.  
2. **DNA & RNA Extraction** – High-quality nucleic acids prepared.  
3. **Genome Sequencing** – Illumina short reads + Nanopore long reads.  
4. **Hybrid Assembly & QC** – Flye + Pilon, evaluated with QUAST and BUSCO.  
5. **Functional Annotation** – MAKER3, dbCAN3 for CAZymes, Geneious for manual curation.  
6. **Comparative Genomics** – OrthoFinder, MCScanX, MAUVE for orthogroups and synteny.  
7. **Phylogenomics** – Tree reconstruction with single-copy orthologs (ML + Bayesian).  
8. **Transcriptomics** – RNA-Seq (12 libraries, 3 conditions × 4 replicates), DEG analysis with DESeq2.  
9. **Integration** – Multi-omics approach to pinpoint candidate enzymes and pathways.  

---

## Summary of Experimental Workflow

| Step | Category | Tools |
|------|----------|-------|
| 1 | Isolation | PDA + Rhodamine B |
| 2 | DNA/RNA Extraction & QC | CTAB, Qubit, Nanodrop |
| 3 | Sequencing | Illumina + Nanopore |
| 4 | Assembly | Flye + Pilon |
| 5 | QC | QUAST, BUSCO |
| 6 | Annotation | MAKER3, dbCAN3, Geneious |
| 7 | Comparative Genomics | OrthoFinder, MCScanX, MAUVE |
| 8 | Phylogeny | MAFFT, RAxML, MrBayes |
| 9 | Transcriptomics | Salmon, DESeq2 |

---

## Detailed Cost Breakdown

| Personnel              | Activity                               | Estimated Cost (€) | Description                              |
| ---------------------- | -------------------------------------- | ------------------ | ---------------------------------------- |
| Wet lab Postdoc salary |                                        | 45,000             | One-year contract for lab activities     |
|                        | Sampling & Isolation                   | 3,000              | Compost collection & fungal culturing    |
|                        | DNA Extraction & QC                    | 2,000              | CTAB reagents, plasticware, gel materials|
|                        | Illumina PE150 Sequencing              | 30,000             | Short-read library preparation & sequencing |
|                        | Oxford Nanopore Sequencing             | 50,000             | Long-read flow cells & library kits      |
|                        | RNA extraction (12 libraries)          | 1,000              | TRIzol & silica column kits              |
|                        | RNA library preparation (12 libraries) | 1,500              | Illumina TruSeq RNA kits                 |
|                        | RNASeq Illumina (12 libraries)         | 2,500              | PE150 sequencing                         |
| Bioinformatics Postdoc |                                        | 45,000             | One-year contract for analysis & reporting |
|                        | Genome Assembly & Polishing            | 0                  | Flye + Pilon (open-source)               |
|                        | Functional Annotation                  | 0                  | MAKER3, dbCAN3, Geneious (manual curation)|
|                        | Comparative Genomics                   | 0                  | OrthoFinder, MAUVE, MCScanX              |
|                        | Non-open source softwares              | 5,000              | Geneious & other licensed tools          |
|                        | Reports, Dissemination                 | 15,000             | Publications, conferences, outreach     |
| **Total**              |                                        | **200,000**        |                                          |

---

## Expected Key Results

- **Genome assemblies**: High contiguity and high completeness (BUSCO >90% expected).  
- **Annotation**: Broad CAZyme repertoire predicted, including esterases, cutinases, lipases.  
- **Comparative genomics**: Orthogroups shared by degraders, absent in non-degraders, highlight candidate genes.  
- **Transcriptomics**: PLA/PHA expected to upregulate hydrolase families; overlap with comparative genomics strengthens candidates.  
- **Integration**: Functional framework of enzymes and evolutionary context for degradation capacity.  

---
## Fungal Candidates in PLA/PHA Waste (Simulated Preview)

Based on **literature and preliminary reports**, the following fungi are typically enriched in PLA/PHA-containing environments:

- **Ascomycota:**  
  - *Aspergillus* spp. (esterase-rich, known PLA degraders)  
  - *Fusarium* spp. (cutinases and depolymerases)  
  - *Penicillium* spp. (broad CAZyme repertoire)  
  - *Purpureocillium* spp. (documented for polyester degradation, also biocontrol)  

- **Basidiomycota:**  
  - *Phanerochaete chrysosporium* (white-rot fungus, strong lignin/peroxide system)  
  - *Trametes* spp. (laccase-rich, oxidative degradation)  

- **Mucoromycota:**  
  - *Rhizopus* spp. (esterase activity, opportunistic degraders)  

---

## Bioeconomy Impact and EEA Policy Alignment

- Supports **EU Circular Economy** by enabling fungal enzymes for PLA/PHA degradation.  
- Provides molecular resources for **industrial composting optimization**.  
- Contributes to reducing **microplastic accumulation**.  
- Aligns with **EEA recommendations** for sustainable plastic alternatives.

**Reference:**  
[EEA — Biodegradable and compostable plastics](https://www.eea.europa.eu/en/analysis/publications/biodegradable-and-compostable-plastics)

---

## Sequencing Technology Overview

| Technology       | Read Length       | Coverage | Mean Q-Score | Advantages                                      | Primary Applications |
| ---------------- | ----------------- | -------- | ------------ | ----------------------------------------------- | -------------------- |
| Illumina NovaSeq | ~150 bp (paired)  | ~101×    | >30          | High-accuracy short reads; optimal for polishing assemblies | SNP calling, differential expression (DEGs) |
| Oxford Nanopore  | ~10–20 kb         | ~35×     | >30          | Long reads; resolves complex repeats and structural variants | De novo assembly, structural genomics |

**Platforms used:**
1. **Illumina NovaSeq 6000**  
   Product info: [Illumina NovaSeq](https://www.illumina.com/systems/sequencing-platforms/novaseq.html)
2. **Oxford Nanopore GridION**  
   Product info: [Nanopore GridION](https://nanoporetech.com/products/gridion)

---

## References

### Scientific Articles

- Menicucci, A., Iacono, S., Ramos, M., Fiorenzani, C., Peres, N. A., et al. (2025). Can whole genome sequencing resolve taxonomic ambiguities in fungi? The case study of Colletotrichum associated with ferns. Frontiers in Fungal Biology, 6, 1540469.
- Khan et al. (2021). Fungal esterases and cutinases: Key players in synthetic polyester degradation. Applied Microbiology and Biotechnology, 105(11), 4445–4463.
- Prasad, P., Varshney, D., & Adholeya, A. (2015). Whole genome annotation and comparative genomic analyses of bio-control fungus Purpureocillium lilacinum. BMC Genomics, 16(1), 1004.
- Ekanayaka, A. H. et al. (2025). Linking metabolic pathways to plastic-degrading fungi: a comprehensive review. J. Fungal Biol.

### Main Tools Used

#### Assembly
- **Flye** – Long-read genome assembler [https://github.com/fenderglass/Flye](https://github.com/fenderglass/Flye)  
- **Pilon** – Assembly polishing with Illumina reads [https://github.com/broadinstitute/pilon](https://github.com/broadinstitute/pilon)  

#### Quality Control (QC)
- **QUAST** – Assembly quality assessment [http://quast.sourceforge.net/](http://quast.sourceforge.net/)  
- **BUSCO** – Benchmarking Universal Single-Copy Orthologs [https://busco.ezlab.org/](https://busco.ezlab.org/)  

#### Annotation
- **MAKER3** – Genome annotation pipeline [http://www.yandell-lab.org/software/maker.html](http://www.yandell-lab.org/software/maker.html)  
- **dbCAN3** – CAZyme annotation server [https://bcb.unl.edu/dbCAN2/](https://bcb.unl.edu/dbCAN2/)  

#### Comparative Genomics
- **OrthoFinder** – Orthogroup inference [https://github.com/davidemms/OrthoFinder](https://github.com/davidemms/OrthoFinder)  
- **MCScanX** – Synteny/collinearity analysis [https://github.com/wyp1125/MCScanX](https://github.com/wyp1125/MCScanX)  

#### Transcriptomics
- **Salmon** – Transcript quantification [https://salmon.readthedocs.io/](https://salmon.readthedocs.io/)  
- **DESeq2** – Differential expression analysis (R/Bioconductor) [https://bioconductor.org/packages/release/bioc/html/DESeq2.html](https://bioconductor.org/packages/release/bioc/html/DESeq2.html)  

### Visualization
- **R / Python** – For plotting and downstream analysis [https://www.r-project.org/](https://www.r-project.org/) | [https://www.python.org/](https://www.python.org/)  
- **Geneious** – Manual curation and visualization [https://www.geneious.com/](https://www.geneious.com/)  

---

## Contact

**Martina Castellucci**  
Email: [martina.castellucci@studio.unibo.it](mailto:martina.castellucci@studio.unibo.it)  
University of Bologna — Applied Genomics Project
