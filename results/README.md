# Multi-omics Analysis Results for *Purpureocillium lilacinum* PLA-C1

This repository contains the main outputs from genome assembly, transcriptome analysis, functional annotation, and comparative genomics.

## 01_BUSCO
- **Assembled genome** (`GCF_023168085.1_PurlilCBS_1.0_genomic.zip`)
- **BUSCO report** (`Galaxy14-[Busco_short_summary].txt`)

## 02_CAZy_annotation
- **Predicted proteome** (`protein.faa`)
- **CAZy annotations** (`overview.txt`, `CAZyme.pep`)

## 03_Orthology_analysis
- **Orthology results** (`Orthogroups.txt.gz`, `Orthologs.txt.gz`)
- **Alignments** (`Multiple_sequence_alignment.fasta.gz`, `single_copy.tar.gz`)

## 04_Phylogeny
- **Phylogenetic tree** (`Species_phylogenetic_tree.nwk.gz`)
- **Tree figure** (PNG generated from the `.nwk` file)

## 05_Transcriptome
- **Transcript nucleotide sequences** (`rna.fna`)
- **Proteins from transcripts** (`gffread_pep.fasta`)

---

## Key figures
- **Venn / UpSet diagrams** → generated from `Orthogroups.txt.gz`
- **Phylogeny** → generated from `Species_phylogenetic_tree.nwk.gz`


