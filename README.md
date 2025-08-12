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
- [Bioeconomy Impact and EEA Policy Alignment](#bioeconomy-impact-and-EEA-policy-alignment)
- [Sequencing Technology Overview](#sequencing-technology-overview)
- [Data Management](#data-management)
- [References](#references)
- [Contact](#contact)

---

## Project Overview

*Purpureocillium lilacinum* strain PLA-C1 was isolated from compost material enriched with biodegradable polylactic acid (PLA) and polyhydroxyalkanoates (PHA). This filamentous fungus exhibits active growth directly on PLA fragments and produces extracellular esterases detectable through rhodamine B plate assays, indicating a high potential for **bioplastic degradation**.

The Applied Genomics study focused on:

1. Generating a **high-quality hybrid de novo genome assembly** from Oxford Nanopore and Illumina sequencing reads obtained directly from the PLA/PHA compost isolate.
2. Performing **comprehensive functional annotation** to identify enzymes involved in the breakdown of synthetic and bio-based polyesters.
3. Conducting **RNA-Seq differential expression analysis** under PLA, PHA, and control growth conditions to determine gene regulation patterns related to polymer degradation.
4. Performing **comparative genomics** and **phylogenomic analyses** with related filamentous fungi to assess unique gene content, synteny, and evolutionary relationships.

---

## Habitat & Ecology

*Purpureocillium lilacinum* is a filamentous ascomycete in the family **Ophiocordycipitaceae**, occurring in a wide range of terrestrial and aquatic habitats, including:

- Agricultural soils and greenhouse substrates  
- Forest and grassland soils  
- Estuarine and river sediments  
- Moist indoor environments  

It is a **multitrophic organism**, functioning as a saprophyte, nematophagous parasite, and endophyte. The species is widely recognized as a **biocontrol agent** against plant-parasitic nematodes and has also been documented in roles related to **bioremediation**.

**Environmental tolerance:**
- Temperature range: 8 °C to ~38 °C
- pH tolerance: 4 to 9
- Resistance to some environmental stresses and chemical treatments

The ability to secrete extracellular hydrolytic enzymes, including cutinases, esterases, and lipases, positions this species as an ideal candidate for **bioplastic waste degradation** in composting environments.

---

## Background and Rationale

Biodegradable polymers such as PLA and PHA are promoted as eco-friendly alternatives to petroleum-based plastics but often degrade incompletely under standard composting conditions. Filamentous fungi represent promising biological agents for bioplastic degradation because they:

- Secrete **extracellular hydrolases** capable of attacking ester bonds in polymers  
- Thrive in **nutrient-limited environments** like compost piles  
- Adapt to diverse carbon sources and environmental stresses  

Despite these advantages, there is a lack of **high-quality genome assemblies** and **integrated transcriptomic datasets** for fungi that can degrade PLA and PHA. This gap hinders the discovery and biotechnological application of relevant enzymes.

---

## Objectives

- Isolate *P. lilacinum* PLA-C1 from compost enriched with PLA and PHA.  
- Produce a **hybrid genome assembly** integrating Nanopore long reads and Illumina short reads from the isolate’s genomic DNA.  
- Annotate genes and predict **polymer-degrading enzymes** through domain-based and homology-based pipelines.  
- Perform **RNA-Seq differential expression analysis** under PLA, PHA, and control conditions to identify genes induced by bioplastic substrates.  
- Conduct **comparative genomics** with related fungal species to reveal lineage-specific genes and syntenic regions.  
- Construct **phylogenomic trees** using single-copy orthologs to clarify evolutionary relationships.

---

## Methodological Pipeline

<img width="515" height="265" alt="image" src="https://github.com/user-attachments/assets/162fb104-93be-4885-978c-9c7e1915781b" />

1. **Sample Collection & Isolation**  
   - Compost originating from PLA/PHA-enriched waste was plated on PDA medium supplemented with 0.5% PLA and Rhodamine B.  
   - Extracellular esterase activity confirmed via rhodamine fluorescence assay under UV illumination.  
   - Identity of the isolate (*Purpureocillium lilacinum* strain PLA-C1) confirmed by ITS rDNA sequencing and BLAST-based phylogenetic placement.  

2. **DNA Extraction & Quality Control**  
   - High-molecular-weight genomic DNA extracted using the CTAB protocol optimized for filamentous fungi.  
   - QC included Qubit fluorometry, Nanodrop spectrophotometry (A260/280, A260/230), and agarose gel electrophoresis (>20 kb bands).  

3. **Genome Sequencing**  
   - Illumina NovaSeq PE150 (~101× coverage) for short reads.  
   - Oxford Nanopore GridION (~35× coverage) for long reads.  

4. **Genome Assembly & Polishing**  
   - Flye v2.9 (Nanopore) for primary assembly.  
   - Pilon v1.24 (3 rounds) for polishing using Illumina reads.  
   - Assembly quality assessed via QUAST (contiguity) and BUSCO (fungi_odb10 lineage dataset).  

5. **Functional Genome Annotation**  
   - Annotation performed with MAKER3, integrating **Augustus** and **GeneMark-ES** predictions with RNA-Seq alignments.  
   - CAZyme prediction via dbCAN3 (HMMER search against CAZy HMM database).  
   - Manual validation and visualization of gene models and annotations performed in **Geneious**.  

6. **Comparative Genomics**  
   - OrthoFinder for orthogroup inference and unique gene identification.  
   - MCScanX and MAUVE for synteny and collinearity analysis.  

7. **Phylogenetic Reconstruction**  
   - 338 single-copy orthologs aligned using MAFFT, conserved regions selected with Gblocks, concatenated using AMAS.  
   - Tree inference with RAxML (maximum likelihood), MrBayes (Bayesian inference), MEGA11 (visualization).  
   - Final trees curated and annotated in **Geneious**.  

8. **Transcriptomic Analysis**  
   - RNA extracted from 3 growth conditions (Control, PLA, PHA) × 4 replicates.  
   - Quantification with Salmon (bias correction enabled) → aggregation with tximport.  
   - Differential expression with DESeq2 (FDR ≤ 0.05, |log2FC| ≥ 2).  

---

## Summary Table of Experimental Workflow

| Step | Category | Tool(s) | Output Files |
|------|----------|---------|--------------|
| 1 | Isolation | PDA+Rhodamine B | – |
| 2 | DNA Extraction/QC | CTAB, Qubit, Nanodrop | – |
| 3 | Sequencing | NovaSeq, GridION | – |
| 4 | Assembly | Flye + Pilon | `assembly.fasta` |
| 5 | QC | QUAST, BUSCO | `Galaxy14-[Busco_summary].txt` |
| 6 | Annotation | MAKER3, GeneMark-ES, Augustus, dbCAN3, Geneious | `genomic.gff`, `protein.faa`, `CAZyme.pep`, `overview.txt` |
| 7 | Comparative Genomics | OrthoFinder, MCScanX, MAUVE | `Orthogroups.txt`, `Orthologs.txt` |
| 8 | Phylogeny | MAFFT, Gblocks, AMAS, RAxML, MrBayes, MEGA11, Geneious | `f1a8dba68a6f43b189057b437429746f-fasta-tree.png` |
| 9 | Transcriptomics | Salmon, tximport, DESeq2 | `rna_counts.tsv`, `deseq2_results.csv` |

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

## Key Results

### Genome Assembly
The hybrid assembly of *P. lilacinum* PLA-C1, generated from Oxford Nanopore and Illumina NovaSeq reads, resulted in a **high-contiguity genome**:

| Metric                 | Value |
|------------------------|-------|
| Genome size            | 38.6 Mb |
| Total length (ungapped)| 38.6 Mb |
| Number of contigs      | 10 |
| N50                    | 5.3 Mb |
| GC content             | 58.5% |
| Genome coverage        | 193× |
| Assembly level         | Contig |
| BUSCO completeness     | **76.3%** (Single-copy: 75.7%, Duplicated: 0.6%, Fragmented: 1.5%, Missing: 22.2%, n=3546) |

**Interpretation:**  
Although the assembly shows excellent contiguity (N50 > 5 Mb, only 10 contigs), the BUSCO completeness (76.3%) suggests the genome may be partially incomplete or that certain lineage-specific genes are absent, possibly due to ecological adaptation to PLA/PHA-enriched environments.


### Functional Annotation
Functional annotation identified **272 CAZymes** (Carbohydrate-Active enZymes) spanning GH, GT, CE, AA, and CBM families. Several hydrolases are overexpressed under PLA conditions.

- **Candidate degradative enzymes:**
  - **Cutinases:** 9 total (5 PLA-induced)
  - **Esterases:** 45 total (12 PLA-induced)
  - **Lipases:** 27 total (7 PLA-induced)
  - **PHA depolymerase-like:** 4 total (2 PLA-induced)

Additional findings:
- **Biosynthetic Gene Clusters (BGCs):** 14 total, including 3 clusters not found in closely related taxa, potentially linked to niche adaptation.
- **AntiSMASH results** revealed secondary metabolite clusters including NRPS-like and terpene synthases, which may contribute to environmental resilience.


### Comparative Genomics (OrthoFinder)
Using OrthoFinder, *P. lilacinum* PLA-C1 was compared to *Penicillium chrysogenum* and *Fusarium solani* proteomes.

**Orthology statistics:**
- Shared orthogroups: **10,312**
- Strain-specific genes: **314**  
- Single-copy clusters: **5,133**
- Percentage of singletons: **21.66%**

**Observation:**  
The high number of unique genes suggests potential novel enzymatic capabilities, possibly linked to PLA/PHA degradation.


### Transcriptomics (DESeq2)
RNA-Seq analysis under three growth conditions (**Control**, **PLA**, **PHA**) identified differentially expressed genes (DEGs) with FDR ≤ 0.05 and |log2FC| ≥ 2:

| Condition comparison | DEGs identified |
|----------------------|-----------------|
| PLA vs Control       | 84 |
| PHA vs Control       | 51 |
| Shared PLA/PHA       | 29 |

**Notable trends:**
- PLA-induced DEGs include multiple cutinases, esterases, and lipases.
- Several DEGs map to **CAZy families GH and CE**, consistent with polymer breakdown.
- PHA response shows fewer DEGs but includes putative PHA depolymerases.


### Summary Insight
The integration of genome assembly, functional annotation, and transcriptomics strongly supports the bioplastic degradation potential of *P. lilacinum* PLA-C1, with **specific hydrolases upregulated under PLA exposure** and unique gene clusters possibly linked to adaptation to synthetic polymer-rich environments.

---

## Bioeconomy Impact and EEA Policy Alignment

The high-quality genome and transcriptome of *Purpureocillium lilacinum* PLA-C1 constitute a **strategic biotechnological resource** for addressing the challenges of bioplastic waste management within the framework of the **EU Circular Economy Action Plan**. This strain, isolated from PLA/PHA-enriched compost, harbors a repertoire of **hydrolase-encoding genes**—including PLA-inducible lipases, cutinases, and esterases—capable of enhancing bioplastic degradation under realistic composting and soil conditions.

These findings directly address the **European Environment Agency (EEA)** concerns regarding biodegradable plastics such as **PLA, PBS, and PBAT**, which often exhibit **slow or incomplete mineralization** in natural environments, contributing to **microplastic accumulation**. By combining **hybrid genome assembly** with **multi-omics analyses**, this study provides:

- **Targeted enzyme discovery** for industrial-scale bioplastic depolymerization.
- **Evidence-based composting optimization**, adjusting parameters to accelerate PLA/PHA breakdown.
- **Mitigation strategies** for microplastic persistence in terrestrial and aquatic ecosystems.
- **Scalable bioprocess integration**, aligning with EEA sustainability and waste management goals.

In doing so, *P. lilacinum* PLA-C1 serves as both a **model organism** for fungal bioplastic degradation and a **practical candidate** for industrial and municipal composting solutions.

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

## Data Management and Accessibility

All raw and processed datasets are systematically organized for reproducibility:

- **`results/`** — BUSCO completeness, CAZy annotation, OrthoFinder orthology results, phylogenetic trees, DESeq2 transcriptomic analysis.
- **`00_Input_data/`** — Assembled genome, predicted proteins, annotated transcripts.

Data are archived in structured directories to facilitate direct integration into downstream pipelines (e.g., enzyme screening, comparative genomics).

---

## References

## References

### Scientific Articles

**Fungal bioplastic degradation and *Purpureocillium lilacinum* studies**

* Prasad, P., Varshney, D., & Adholeya, A. (2015). Whole genome annotation and comparative genomic analyses of bio-control fungus *Purpureocillium lilacinum*. **BMC Genomics**, 16:1004. DOI: [10.1186/s12864-015-2229-2](https://doi.org/10.1186/s12864-015-2229-2)  
* Xie, J.-L. et al. (2016). Genome and transcriptome sequences reveal the specific parasitism of the nematophagous *Purpureocillium lilacinum* 36-1. **Frontiers in Microbiology**, 7:1084. DOI: [10.3389/fmicb.2016.01084](https://doi.org/10.3389/fmicb.2016.01084)  
* Li, Y. et al. (2024). Revealing the metabolic potential and environmental adaptation of piezotolerant *Purpureocillium lilacinum* FDZ8Y1 from the Mariana Trench. **Front. Microbiol.** DOI: [10.3389/fmicb.2024.1474180](https://www.frontiersin.org/journals/microbiology/articles/10.3389/fmicb.2024.1474180/full)  
* Tseng, W. S. et al. (2023). Poly(butylene adipate-co-terephthalate) biodegradation by *Purpureocillium lilacinum* strain BA1S. **Appl Microbiol Biotechnol**. DOI: [10.1007/s00253-023-12566-9](https://doi.org/10.1007/s00253-023-12566-9)  
* Urbanek, A.K. et al. (2018). Biodegradation of plastics by fungal communities – opportunities and limitations. **Appl Microbiol Biotechnol**. DOI: [10.1007/s00253-018-9271-y](https://doi.org/10.1007/s00253-018-9271-y)  
* Harms, H. et al. (2021). Plastics in the environment – fungal enzymes to the rescue? **Biotechnol Adv**. DOI: [10.1016/j.biotechadv.2021.107712](https://doi.org/10.1016/j.biotechadv.2021.107712)  
* Ekanayaka, A. H. et al. (2025). Linking metabolic pathways to plastic-degrading fungi: a comprehensive review. **J. Fungal Biol.** DOI: [10.3390/jof11050378](https://www.mdpi.com/2309-608X/11/5/378)  
* Lo Giudice, A. et al. (2023). Purpureocillium lilacinum for biocontrol and bioremediation: insights from environmental and applied microbiology. **Processes**, 11(10):3445. DOI: [10.3390/pr11103445](https://doi.org/10.3390/pr11103445)  

### Bioinformatics Tools & Methods

**Genome Assembly & Polishing**  
* Flye – Long-read genome assembler.  
* Pilon – Assembly polishing with short reads.  
* QUAST – Assembly quality assessment.  
* BUSCO – Genome completeness evaluation.  

**Genome Annotation & Functional Genomics**  
* MAKER3 – Genome annotation pipeline integrating ab initio predictors and transcript evidence.  
* Augustus – Ab initio gene prediction.  
* GeneMark-ES – Self-training gene prediction.  
* dbCAN3 – Carbohydrate-active enzyme annotation.  
* Geneious – Manual curation and secondary validation of annotations.  

**Comparative Genomics & Phylogenetics**  
* OrthoFinder – Orthogroup inference.  
* MCScanX – Synteny and collinearity detection.  
* MAUVE – Whole-genome alignment.  
* MAFFT – Multiple sequence alignment.  
* Gblocks – Selection of conserved sequence blocks.  
* AMAS – Alignment concatenation.  
* RAxML – Maximum likelihood phylogenetic analysis.  
* MrBayes – Bayesian phylogenetic inference.  
* MEGA11 – Tree visualization and evolutionary analysis.  
* Geneious – Tree visualization and annotation.  

**Transcriptomics & DEG Analysis**  
* Salmon – Transcript quantification.  
* tximport – Aggregation of transcript counts.  
* DESeq2 – Differential expression analysis.

### Data Sources & Platforms
* Galaxy Europe – [https://usegalaxy.eu/](https://usegalaxy.eu/) (workflow execution and analysis environment)  
* CAZy Database – [http://www.cazy.org/](http://www.cazy.org/) (reference for carbohydrate-active enzymes)  
* KEGG – [https://www.genome.jp/kegg/](https://www.genome.jp/kegg/) (pathway annotation)  

---

## Contact
Martina Castellucci  
Email: [martina.castellucci@studio.unibo.it](mailto:martina.castellucci@studio.unibo.it)  
University of Bologna — Applied Genomics Project
