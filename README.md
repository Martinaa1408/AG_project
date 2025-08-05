# Hybrid Genome Assembly and Integrative Multi-Omics of *Purpureocillium lilacinum* PLA-C1 Reveals Enzymatic Potential for Bioplastic Degradation

[![GitHub](https://img.shields.io/badge/GitHub-Repository-181717?logo=github&logoColor=white)](https://github.com/YourUser/PLA-C1-GenomeProject)
[![NCBI BioProject](https://img.shields.io/badge/NCBI-BioProject-336699?logo=databricks&logoColor=white)](https://www.ncbi.nlm.nih.gov/bioproject/)
[![ENA Sequencing](https://img.shields.io/badge/ENA-Reads%20%26%20RNASeq-FF6600?logo=databricks&logoColor=white)](https://www.ebi.ac.uk/ena/browser/home)
[![Illumina NovaSeq](https://img.shields.io/badge/Illumina-NovaSeq%20PE150-FF0000?logo=illumina&logoColor=white)](https://www.illumina.com/systems/sequencing-platforms/novaseq.html)
[![Oxford Nanopore](https://img.shields.io/badge/Nanopore-GridION-0078D4?logo=nanopore&logoColor=white)](https://nanoporetech.com/products/gridion)
[![Contact](https://img.shields.io/badge/Contact-Email-EA4335?logo=gmail&logoColor=white)](mailto:martina.castellucci@studio.unibo.it)
[![PDF Report](https://img.shields.io/badge/Download-Final_Report-lightgrey?logo=adobeacrobatreader&logoColor=red)](link_al_pdf)


## Table of Contents

* [Project Overview](#project-overview)
* [Background and Rationale](#background-and-rationale)
* [Objectives](#objectives)
* [Methodological Pipeline](#methodological-pipeline)
* [Summary Table of Experimental Workflow](#summary-table-of-experimental-workflow)
* [Detailed Cost Breakdown](#detailed-cost-breakdown)
* [Expected Results and Bioeconomy Impact](#expected-results-and-bioeconomy-impact)
* [Sequencing Technology Overview](#sequencing-technology-overview)
* [Data Management](#data-management)
* [References](#references)
* [Contact](#contact)

---

## Project Overview

*Purpureocillium lilacinum* strain PLA-C1 was isolated from compost enriched with biodegradable polylactic acid (PLA). This environmental filamentous fungus shows growth directly on PLA fragments and produces extracellular esterases detectable by rhodamine B screening, indicating a potential for **bioplastic degradation**.

This **Applied Genomics project** focuses on:

1. **Hybrid de novo genome assembly** (~41.6 Mb, N50 1.32 Mb, 119 contigs, 98.3% BUSCO) using **Oxford Nanopore** and **Illumina** data.  
2. **Functional annotation** to identify **candidate degradative genes** (cutinases, esterases, lipases, PHA depolymerases) and **14 biosynthetic gene clusters (BGCs)**, 3 unique to this strain.  
3. **Transcriptomic profiling (RNA-Seq)** under PLA, PHA, and control conditions to detect **differentially expressed genes** involved in polymer degradation.  
4. **Comparative and phylogenomic analyses** to determine the environmental adaptation and unique gene content relevant to bioplastic degradation.

This genomic resource provides the **foundation for functional enzyme discovery** and supports **bioeconomic applications** in **compost-based bioplastic waste management**.

---

## Background and Rationale

Plastics, including “biodegradable” ones such as PLA and PHA, often persist under real composting conditions, creating environmental concerns.  

Filamentous fungi are promising candidates for **bioplastic biodegradation** because they:

- Secrete **extracellular hydrolases** (esterases, cutinases, lipases)  
- Can adapt to **nutrient-limited environments** like compost  
- Exhibit **metabolic flexibility** to degrade complex polymers

**Gap in knowledge:**  
Few **high-quality genome assemblies** and **transcriptomic datasets** exist for environmental fungi capable of PLA degradation, limiting enzyme discovery for applied biotechnology.

**Aim of this project:**  
Combine **hybrid genomics and transcriptomics** to **characterize the enzymatic potential** of *P. lilacinum* PLA-C1 for **applied bioplastic degradation**.

---

## Objectives

* Isolate and characterize **P. lilacinum PLA-C1** from compost enriched with PLA.
* Produce a **high-quality hybrid de novo genome assembly** and evaluate completeness (BUSCO).  
* Perform **functional annotation** to identify candidate degradative enzymes and BGCs.  
* Profile **gene expression via RNA-Seq** under PLA/PHA vs control to detect **induced genes**.  
* Compare the genome with related fungi to detect **lineage-specific genes and synteny**.  
* Generate **genomic and transcriptomic resources** for **applied environmental biotechnology**.

---

## Methodological Pipeline

<img width="515" height="265" alt="image" src="https://github.com/user-attachments/assets/162fb104-93be-4885-978c-9c7e1915781b" />


### 1. Sample Collection and Fungal Isolation

* Compost samples with partially degraded PLA cups (Bologna, Italy).  
* Isolation on PDA + 0.5% PLA + 0.02% Rhodamine B (esterase screen).  
* ITS rDNA sequencing confirmed **Purpureocillium lilacinum PLA-C1**.

### 2. DNA Extraction and QC

* CTAB-based protocol for filamentous fungi.  
* QC: Qubit (130 ng/µL), Nanodrop (A260/280 1.85), agarose gel (band >20 kb).  

### 3. Genome Sequencing

* **Illumina NovaSeq PE150 (~100X)** → high accuracy, polishing, RNA-Seq  
* **Oxford Nanopore GridION (~30–50X)** → long reads for assembly  

### 4. Assembly and Quality Assessment

* Flye (long-read assembly) + Pilon (short-read polishing, 3 rounds)  
* QC via QUAST and BUSCO (fungi_odb10)  

**Expected Assembly Metrics**:

| Feature             | Value      |
|---------------------|-----------|
| Genome size         | ~41.6 Mb  |
| Contigs             | 119       |
| N50                 | 1.32 Mb   |
| GC content          | 50.2%     |
| BUSCO completeness  | 98.3%     |

### 5. Transcriptome Analysis

* 3 conditions × 4 replicates (12 RNA-Seq libraries):  
  - Control (MM + glucose 1%)  
  - PLA (MM + 1% PLA fragments)  
  - PHA (MM + 1% PHB)  
* RNA integrity RIN ≥ 8.5, ~20M PE150 reads/sample  
* Pipeline: **STAR → featureCounts → DESeq2** (log2FC>2, FDR<0.05)

**Differential Expression**:  
- 84 DEGs in PLA  
- 29 DEGs shared PLA/PHA

### 6. Genome Annotation

* **MAKER3** + Augustus + GeneMark + RNA-Seq evidence  
* Functional annotation: InterProScan, Pfam, eggNOG, KEGG, CAZy, AntiSMASH  

**Candidate Enzymes**:

| Enzyme Family         | Domain        | Genes | PLA DEGs |
|-----------------------|---------------|-------|---------|
| Cutinase              | IPR000675      | 9     | 5       |
| Esterase              | IPR000379      | 45    | 12      |
| Lipase                | IPR000734      | 27    | 7       |
| PHA depolymerase-like | UniProt HMM    | 4     | 2       |

### 7. Comparative Genomics

* Reference: *Talaromyces purpureogenus*, *Paecilomyces variotii*, *Fusarium solani*  
* OrthoFinder → 10,312 shared orthogroups, 314 unique genes  
* MAUVE & MCScanX → synteny and unique segments

### 8. Phylogenetic Analysis

* 338 single-copy orthologs (MAFFT → Gblocks → AMAS)  
* Trees: RAxML + MrBayes + MEGA11 → PLA-C1 forms distinct clade

---

## Summary Table of Experimental Workflow

| Step | Category              | Tool/Protocol                         | Output                                 |
| ---- | --------------------- | ------------------------------------- | -------------------------------------- |
| 1    | Sample isolation      | PLA compost plating (PDA+Rhodamine B) | Fungal isolate (PLA-C1)                |
| 2    | DNA extraction        | CTAB protocol                         | High-MW DNA (>20 kb)                   |
| 3    | DNA QC                | Nanodrop, Qubit, gel electrophoresis  | DNA purity & concentration             |
| 4    | Illumina sequencing   | NovaSeq PE150                         | ~100X PE reads                         |
| 5    | Nanopore sequencing   | GridION (LSK-109)                     | Long reads (~30–50X)                   |
| 6    | Assembly              | Flye + Pilon                          | `assembly.fasta`, `quast_report.txt`   |
| 7    | RNA-Seq               | Illumina PE150                        | 12 libraries                           |
| 8    | Expression analysis   | STAR + DESeq2                         | `rna_counts.tsv`, `deseq2_results.csv` |
| 9    | Gene prediction       | MAKER3 + Augustus/GeneMark            | `annotation.gff3`, `transcripts.fasta` |
| 10   | Functional annotation | InterPro, Pfam, KEGG, AntiSMASH       | Annotated domains, enzymes, BGCs       |
| 11   | Comparative genomics  | OrthoFinder + MAUVE/MCScanX           | `orthogroups.tsv`, synteny plots       |
| 12   | Phylogenetics         | MAFFT, RAxML, MrBayes, MEGA11         | `phylogenetic_tree.nwk`                |

---

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

## Expected Results and Bioeconomy Impact

* High-quality **annotated genome and transcriptome** of *P. lilacinum* PLA-C1  
* **Candidate enzyme catalog** for PLA/PHA degradation  
* **Unique BGCs and lineage-specific genes** revealed by comparative genomics  
* **Differential expression data** confirming enzyme induction under PLA exposure  
* Genomic resource supporting **bioplastic composting and circular bioeconomy**

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

* Raw reads → **ENA BioProject (TBA)**  
* Assembly + annotations + expression matrices → **GitHub repository**  
* **Snakemake workflows** and **Jupyter notebooks** for reproducibility  
* **MIxS-compliant metadata** for FAIR data sharing

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
