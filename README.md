# Hybrid Genome Assembly and Integrative Multi-Omics of *Purpureocillium lilacinum* PLA-C1 Reveals Enzymatic Potential for Bioplastic Degradation

[![Fungal Isolate](https://img.shields.io/badge/Fungal_Isolate-PLA--C1-6A0DAD?logo=leaflet&logoColor=white)]()
[![DNA Extraction](https://img.shields.io/badge/DNA-CTAB_Protocol-228B22?logo=dna&logoColor=white)]()
[![QC](https://img.shields.io/badge/QC-Nanodrop%2C_Qubit%2C_Gel-1E90FF?logo=zoom&logoColor=white)]()
[![Hybrid Assembly](https://img.shields.io/badge/Assembly-Flye+Pilon-FFD700?logo=files&logoColor=white)]()
[![RNA-Seq](https://img.shields.io/badge/RNA--Seq-Illumina_12_Libraries-FF69B4?logo=rna&logoColor=white)]()
[![Expression Analysis](https://img.shields.io/badge/Expression-STAR+DESeq2-8A2BE2?logo=graph&logoColor=white)]()
[![Gene Prediction](https://img.shields.io/badge/Gene_Prediction-MAKER3+Augustus-20B2AA?logo=dna&logoColor=white)]()
[![Functional Annotation](https://img.shields.io/badge/Functional_Annotation-InterPro%2C_KEGG%2C_AntiSMASH-32CD32?logo=network&logoColor=white)]()
[![Comparative Genomics](https://img.shields.io/badge/Comparative_Genomics-OrthoFinder%2C_MCScanX-DA70D6?logo=git-branch&logoColor=white)]()
[![Phylogenetics](https://img.shields.io/badge/Phylogenetics-MAFFT%2C_RAxML%2C_MEGA11-DC143C?logo=tree&logoColor=white)]()

---

## Table of Contents
- [Project Overview](#project-overview)
- [Habitat & Ecology](#habitat--ecology)
- [Background and Rationale](#background-and-rationale)
- [Objectives](#objectives)
- [Methodological Pipeline](#methodological-pipeline)
- [Summary Table of Experimental Workflow](#summary-table-of-experimental-workflow)
- [Detailed Cost Breakdown](#detailed-cost-breakdown)
- [Key Results](#key-results)
- [Expected Bioeconomy Impact](#expected-bioeconomy-impact)
- [Sequencing Technology Overview](#sequencing-technology-overview)
- [Data Management](#data-management)
- [References](#references)
- [Contact](#contact)

---

## Project Overview

*Purpureocillium lilacinum* strain PLA-C1 was isolated from compost enriched with biodegradable polylactic acid (PLA) and polyhydroxyalkanoates (PHA). This filamentous fungus demonstrates growth directly on PLA fragments and produces extracellular esterases detectable by rhodamine B screening, indicating **bioplastic degradation potential**.

The **Applied Genomics project** aimed to:
1. Generate a **hybrid de novo genome assembly** (Nanopore + Illumina).
2. Perform **functional annotation** to identify polymer-degrading enzymes.
3. Conduct **RNA-Seq** differential expression analysis under PLA, PHA, and control conditions.
4. Perform **comparative genomics and phylogeny** with related fungal species.

---

## Habitat & Ecology

*Purpureocillium lilacinum* is a filamentous fungus in the **Ophiocordycipitaceae** family, widely distributed in soils, sediments, rhizospheres, and decaying organic matter. It occurs in:
- Agricultural soils
- Forest and grassland ecosystems
- Estuarine sediments
- Moist indoor environments

It is **multitrophic** — acting as a saprophyte, nematophagous parasite, or endophyte — and is used as a **biocontrol agent** against plant-pathogenic nematodes.

**Environmental tolerance**:
- Temperature: 8 °C to ~38 °C
- pH: 4 to 9
- Tolerant to various environmental stresses and some disinfectants

Its enzyme secretion profile (cutinases, esterases, lipases) makes it an ideal candidate for **bioremediation** and **bioplastic degradation** in compost-based waste management systems.

---

## Background and Rationale

Although PLA and PHA are marketed as biodegradable, they degrade slowly under natural composting conditions. Fungi like *P. lilacinum*:
- Secrete **extracellular hydrolases**
- Tolerate harsh environments
- Can metabolize complex carbon sources

This study fills the **gap in high-quality genomic and transcriptomic resources** for fungi capable of bioplastic degradation.

---

## Objectives

- Isolate *P. lilacinum* PLA-C1 from PLA/PHA-enriched compost.
- Generate a **high-quality hybrid genome assembly**.
- Identify **candidate degradative enzymes** via functional annotation.
- Conduct **RNA-Seq differential expression** under PLA and PHA.
- Compare genome structure with related taxa.
- Build phylogenomic trees from **single-copy orthologs**.

---

## Methodological Pipeline

<img width="515" height="265" alt="image" src="https://github.com/user-attachments/assets/162fb104-93be-4885-978c-9c7e1915781b" />

1. **Sample Collection & Isolation**  
   - PDA + 0.5% PLA + Rhodamine B  
   - ITS sequencing confirmation

2. **DNA Extraction & QC**  
   - CTAB protocol  
   - Qubit, Nanodrop, agarose gel

3. **Genome Sequencing**  
   - Illumina NovaSeq PE150 (~100X coverage)  
   - Oxford Nanopore GridION (~35X coverage)

4. **Assembly & Polishing**  
   - Flye → Pilon (3 rounds)  
   - QUAST, BUSCO (fungi_odb10) for QC

5. **Functional Annotation**  
   - MAKER3 + Augustus + GeneMark  
   - InterProScan, Pfam, eggNOG, KEGG, dbCAN (CAZy), AntiSMASH

6. **Comparative Genomics**  
   - OrthoFinder, MCScanX, MAUVE

7. **Phylogeny**  
   - MAFFT → Gblocks → AMAS  
   - Tree inference: RAxML, MrBayes, MEGA11

8. **Transcriptomics**  
   - RNA extraction from 3 conditions × 4 replicates  
   - Salmon quantification → tximport → DESeq2

---

## Summary Table of Experimental Workflow

| Step | Category | Tool(s) | Output Files |
|------|----------|---------|--------------|
| 1 | Isolation | PDA+Rhodamine B | – |
| 2 | DNA Extraction/QC | CTAB, Qubit, Nanodrop | – |
| 3 | Sequencing | NovaSeq, GridION | `illumina_reads.fastq.gz`, `nanopore_reads.fastq.gz` |
| 4 | Assembly | Flye + Pilon | `assembly.fasta` |
| 5 | QC | QUAST, BUSCO | `busco_summary.txt` |
| 6 | Annotation | MAKER3, InterPro, AntiSMASH, dbCAN | `genomic.gff`, `protein.faa`, `CAZyme.pep` |
| 7 | Comparative Genomics | OrthoFinder, MCScanX | `Orthogroups.txt`, `Orthologs.txt` |
| 8 | Phylogeny | MAFFT, RAxML, MrBayes | `phylogenetic_tree.png` |
| 9 | Transcriptomics | Salmon, DESeq2 | `rna_counts.tsv`, `deseq2_results.csv` |

---

## Detailed Cost Breakdown

| Category | Cost (€) | Notes |
|----------|----------|-------|
| Personnel (wet lab, bioinformatics) | 90,000 | 1-year |
| Sampling & Isolation | 3,000 | Compost handling |
| DNA extraction & QC | 2,000 | Reagents |
| Illumina sequencing | 30,000 | Short reads |
| Nanopore sequencing | 50,000 | Long reads |
| RNA-Seq prep & sequencing | 5,000 | 12 libraries |
| Software licenses | 5,000 | Non-open-source |
| Dissemination | 15,000 | Publications, conferences |
| **Total** | **200,000** | |

---

## Key Results

### Genome Assembly
- Genome size: **38.6 Mb**
- Contigs: **10**
- N50: **5.3 Mb**
- GC content: **58.5%**
- BUSCO completeness: **76.3%**

### Functional Annotation
- **Candidate degradative enzymes**:
  - Cutinases: 9 genes (5 PLA-induced)
  - Esterases: 45 genes (12 PLA-induced)
  - Lipases: 27 genes (7 PLA-induced)
  - PHA depolymerase-like: 4 genes (2 PLA-induced)
- **Biosynthetic Gene Clusters (BGCs)**: 14 total, 3 unique

### Comparative Genomics (OrthoFinder)
- 10,312 shared orthogroups  
- 314 strain-specific genes

### Transcriptomics (DESeq2)
- PLA: 84 DEGs  
- PHA: 65 DEGs  
- PLA/PHA shared: 29 DEGs

---

## Expected Bioeconomy Impact

This genomic and transcriptomic resource can:
- Guide **enzyme discovery** for industrial bioplastic degradation
- Inform **composting optimization** for PLA/PHA waste
- Support the **EU Circular Economy Action Plan**

---

## Sequencing Technology Overview

| Technology       | Read Length       | Coverage | Phred Score | Advantages                            | Typical Use          |
| ---------------- | ----------------- | -------- | ----------- | ------------------------------------- | -------------------- |
| Illumina NovaSeq | ~150 bp (paired)  | ~100X    | >30         | High accuracy, polishing assemblies   | SNPs, DEGs           |
| Oxford Nanopore  | ~10–20 kb         | ~30–50X  | >30         | Long reads for complex regions        | De novo assembly     |

1. **Illumina NovaSeq 6000**
   - Product page: [https://www.illumina.com/systems/sequencing-platforms/novaseq.html](https://www.illumina.com/systems/sequencing-platforms/novaseq.html)

2. **Oxford Nanopore GridION**
   - Product page: [https://nanoporetech.com/products/gridion](https://nanoporetech.com/products/gridion)

---

## Data Management

All input, intermediate, and output files are stored in:
- `results/` — BUSCO, CAZy, OrthoFinder, phylogeny, transcriptomics
- `00_Input_data/` — genome, proteins, transcripts

---

## References



---

## Contact
Martina Castellucci  
Email: [martina.castellucci@studio.unibo.it](mailto:martina.castellucci@studio.unibo.it)  
University of Bologna — Applied Genomics Project





## Detailed Cost Breakdown

| Personnel              | Activity                               | Estimated Cost (€) | Description                              |
| ---------------------- | -------------------------------------- | ------------------ | ---------------------------------------- |
| Wet lab Postdoc salary |                                        | 45,000             | One year project                         |
|                        | Sampling & Isolation                   | 3,000              | Compost handling, fungal plating         |
|                        | DNA Extraction & QC                    | 2,000              | Reagents & consumables                   |
|                        | Illumina PE150 Sequencing              | 30,000             | Short-read library + sequencing          |
|                        | Oxford Nanopore Sequencing             | 50,000             | Long-read sequencing                     |
|                        | RNA extraction (12 libraries)          | 1,000              |                                          |
|                        | RNA library preparation (12 libraries) | 1,500              |                                          |
|                        | RNASeq Illumina (12 libraries)         | 2,500              |                                          |
| Bioinformatics Postdoc |                                        | 45,000             | One year project                         |
|                        | Genome Assembly & Polishing            | 0                  | Flye + Pilon                             |
|                        | Functional Annotation                  | 0                  | MAKER3, InterProScan, AntiSMASH          |
|                        | Comparative Genomics                   | 0                  | OrthoFinder, MAUVE, MCScanX              |
|                        | Non-open source softwares              | 5,000              | Genious & others                         |
|                        | Reports, Dissemination                 | 15,000             | Publications & Conferences               |
| **Total**              |                                        | **200,000**        |                                          |

---



## Environmental Relevance and EEA Alignment

The study aligns with the European Environment Agency (EEA) guidelines on biodegradable and compostable plastics (EEA, 2020).  
According to the EEA, bio-based plastics such as **PBS, PBAT, and PLA** often exhibit **slow or incomplete biodegradation** in soil and aquatic environments, leading to potential microplastic accumulation.  

The **hybrid genome assembly and integrative multi-omics** analysis of *Purpureocillium lilacinum* PLA‑C1 identifies **lipases, cutinases, and esterases** that could **enhance bioplastic degradation** under realistic environmental conditions.  
This connection highlights the potential contribution of such biotechnological studies to the **EU circular economy goals** and to the sustainable management of biodegradable plastic waste.

Reference: [EEA — Biodegradable and compostable plastics](https://www.eea.europa.eu/en/analysis/publications/biodegradable-and-compostable-plastics)

---

## References

### Scientific Articles

**Fungal bioplastic degradation and *Purpureocillium lilacinum* studies**

* Prasad, P., Varshney, D., & Adholeya, A. (2015). Whole genome annotation and comparative genomic analyses of bio-control fungus *Purpureocillium lilacinum*. **BMC Genomics**, 16:1004. DOI: [10.1186/s12864-015-2229-2](https://doi.org/10.1186/s12864-015-2229-2)  
* Xie, J.-L. et al. (2016). Genome and transcriptome sequences reveal the specific parasitism of the nematophagous *Purpureocillium lilacinum* 36‑1. **Frontiers in Microbiology**, 7:1084. DOI: [10.3389/fmicb.2016.01084](https://doi.org/10.3389/fmicb.2016.01084)  
* Li, Y. et al. (2024). Revealing the metabolic potential and environmental adaptation of piezotolerant *Purpureocillium lilacinum* FDZ8Y1 from the Mariana Trench. **Front. Microbiol.** DOI: [10.3389/fmicb.2024.1474180](https://www.frontiersin.org/journals/microbiology/articles/10.3389/fmicb.2024.1474180/full)  
* Tseng, W. S. et al. (2023). Poly(butylene adipate‑co‑terephthalate) biodegradation by *Purpureocillium lilacinum* strain BA1S. **Appl Microbiol Biotechnol**. DOI: [10.1007/s00253-023-12566-9](https://doi.org/10.1007/s00253-023-12566-9)  
* Urbanek, A.K. et al. (2018). Biodegradation of plastics by fungal communities – opportunities and limitations. **Appl Microbiol Biotechnol**. DOI: [10.1007/s00253-018-9271-y](https://doi.org/10.1007/s00253-018-9271-y)  
* Harms, H. et al. (2021). Plastics in the environment – fungal enzymes to the rescue? **Biotechnol Adv**. DOI: [10.1016/j.biotechadv.2021.107712](https://doi.org/10.1016/j.biotechadv.2021.107712)  
* Ekanayaka, A. H. et al. (2025). Linking metabolic pathways to plastic-degrading fungi: a comprehensive review. **J. Fungal Biol.** DOI: [10.3390/jof11050378](https://www.mdpi.com/2309-608X/11/5/378)

### Bioinformatics Tools & Methods

**Genome Assembly & Polishing**
* Flye – Long-read genome assembler. [GitHub](https://github.com/fenderglass/Flye)  
* Pilon – Assembly polishing using Illumina reads. [GitHub](https://github.com/broadinstitute/pilon)  
* QUAST – Quality assessment of genome assemblies. [http://quast.sourceforge.net](http://quast.sourceforge.net)  

**Genome Annotation & Functional Genomics**
* MAKER3 – Genome annotation framework integrating ab initio and RNA-Seq evidence.  
* InterProScan 5 – Genome-scale protein function classification. DOI: [10.1093/bioinformatics/btu031](https://doi.org/10.1093/bioinformatics/btu031)  
* eggNOG-mapper – Orthology and functional annotation. [http://eggnog-mapper.embl.de](http://eggnog-mapper.embl.de)  
* Pfam & CAZy – Protein domains and carbohydrate-active enzymes.  

**Transcriptomics & DEG Analysis**
* STAR – RNA-Seq aligner. DOI: [10.1093/bioinformatics/bts635](https://doi.org/10.1093/bioinformatics/bts635)  
* featureCounts – Read summarization. DOI: [10.1093/bioinformatics/btt656](https://doi.org/10.1093/bioinformatics/btt656)  
* DESeq2 – Differential expression analysis. DOI: [10.1186/s13059-014-0550-8](https://doi.org/10.1186/s13059-014-0550-8)  

**Comparative Genomics & Phylogenetics**
* OrthoFinder – Orthology inference. DOI: [10.1186/s13059-019-1832-y](https://doi.org/10.1186/s13059-019-1832-y)  
* MAUVE – Whole-genome alignment and synteny. [http://darlinglab.org/mauve/mauve.html](http://darlinglab.org/mauve/mauve.html)  
* MCScanX – Detection of synteny and collinearity. [GitHub](https://github.com/wyp1125/MCScanX)  
* RAxML – Phylogenetic analysis (ML). [https://github.com/stamatak/standard-RAxML](https://github.com/stamatak/standard-RAxML)  
* MrBayes – Bayesian phylogenetic inference. [https://nbisweden.github.io/MrBayes/](https://nbisweden.github.io/MrBayes/)  
* MEGA11 – Phylogenetic tree visualization and evolutionary analysis. [https://www.megasoftware.net/](https://www.megasoftware.net/)

### Databases & Taxonomic Resources

* [NCBI Taxonomy – *Purpureocillium lilacinum* (TaxID 123399)](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=123399)  
* [UniProt Taxonomy 123399](https://www.uniprot.org/taxonomy/123399)  
* [MycoBank – *Purpureocillium lilacinum*](https://www.mycobank.org/page/Name%20details/305200)  
* [GBIF Occurrence Data – *Purpureocillium lilacinum*](https://www.gbif.org/species/5241926)  
* [CAZy Database](http://www.cazy.org/) – Carbohydrate-active enzyme families.  
* [KEGG](https://www.genome.jp/kegg/) – Pathways and metabolic annotations.  

---

## Contact

Martina Castellucci  
Email: [martina.castellucci@studio.unibo.it](mailto:martina.castellucci@studio.unibo.it)  
Applied Genomics Project
University of Bologna
