# Applied Genomics Project — simulated workflow pipeline

This README provides a **single, complete, narrative workflow** of the project.  
Each step includes tools, requirements, mock Bash procedures, expected outcomes, and theoretical links.  
All commands are **illustrative only** — they represent the expected pipeline, not real runs.

---

## 0) Prerequisites & Environment

**Concept:** Ensure reproducibility with isolated environments, clear folder layout, and explicit compute resources.

**Tools required:**  
- `mamba/conda`, `git`, `curl/wget`  
- `python >=3.10`, `java >=11`, `R >=4.2`  
- Assemblers and QC: `flye`, `pilon`, `quast`, `busco`, `minimap2`, `bwa`, `samtools`  
- Annotation: `maker`, `augustus`, `genemark-es`, `interproscan`, `hmmer`, `diamond`, `blast`  
- RNA-Seq: `salmon`, `tximport` (R), `deseq2` (R)  
- Comparative: `orthofinder`, `mafft`, `gblocks`, `raxml`, `mrbayes`, `MCScanX`  
- Utilities: `fastqc`, `multiqc`, `trimmomatic`, `fastp`, `seqkit`, `parallel`  

**Suggested folder layout:**
AG_project/
├─ data/ # raw FASTQ, reference, GTF/GFF
├─ qc/ # read QC
├─ assembly/ # Flye + Pilon outputs
├─ annotation/ # MAKER3 / AUGUSTUS / dbCAN3
├─ rna/ # Salmon quants + DESeq2
├─ comparative/ # OrthoFinder / MCScanX / phylogeny
├─ logs/ # HPC logs (SLURM/SGE)
└─ envs/ # conda env files


**Environment creation (example with mamba):**
```bash
mamba create -y -n agproj \
  python=3.11 r-base=4.3 openjdk=17 \
  fastqc trimmomatic fastp pigz multiqc \
  flye quast busco minimap2 samtools bwa pilon \
  maker augustus hmmer diamond blast interproscan \
  orthofinder mafft emboss gblocks raxml mrbayes \
  seqkit parallel
mamba activate agproj
```
Compute needs:

RAM: assembly ≥128 GB; comparative ≥256 GB; RNA-Seq ≤64 GB

CPU: 16–40 cores typical

GPU: optional for Nanopore basecalling (NVIDIA V100/A100)

## 1) Sampling & Screening

Concept: Compost enriched with PLA/PHA is plated. Screening uses PLA + Rhodamine B to identify isolates producing esterases.

Expected: ~20 isolates, ~5 strong degraders.
Theory link: Environmental selection enriches functional degraders.
Galaxy alternative: not applicable (wet-lab step).

# Mock note: colonies with fluorescent halos under UV are retained

## 2) DNA/RNA Extraction & QC

Concept: High molecular weight DNA and intact RNA are critical for NGS.

Quality thresholds:

DNA: A260/280 ~1.8; A260/230 ~2.0; >20 kb fragments

RNA: RIN >8; A260/280 ~2.0

Procedure (mock):

```bash
mkdir -p qc
fastqc -t 8 data/illumina/*fastq.gz data/nanopore/*fastq.gz -o qc/
multiqc qc/ -o qc/
```

Expected: sequencing-grade DNA and RNA.
Theory link: poor input compromises downstream results.
Galaxy alternative: FastQC, MultiQC.

## 3) Hybrid Sequencing

Concept: Illumina NovaSeq provides accurate short reads; ONT GridION generates long reads for contiguity.

Procedure (mock):

```bash
# Basecalling Nanopore reads (GPU, optional)
dorado basecaller dna_r10.4.1_e8.2_400bps_sup@v4.3.0 /path/to/pod5/ > data/nanopore/reads.fastq

# Illumina trimming
fastp -i data/illumina/R1.fastq.gz -I data/illumina/R2.fastq.gz \
      -o data/illumina/R1.trim.fastq.gz -O data/illumina/R2.trim.fastq.gz \
      --detect_adapter_for_pe --thread 16

```
Expected: ~100× Illumina + 30–40× Nanopore coverage.
Theory link: hybrid strategies balance accuracy and contiguity.
Galaxy alternative: Fastp, Cutadapt, NanoPlot.

## 4) Assembly & Polishing

Concept: Flye assembles Nanopore long reads; Pilon polishes with Illumina. QC with QUAST and BUSCO.

Procedure (mock):

```bash
flye --nano-raw data/nanopore/reads.fastq --out-dir assembly/flye --threads 32

bwa index assembly/flye/assembly.fasta
bwa mem -t 32 assembly/flye/assembly.fasta data/illumina/R1.trim.fastq.gz data/illumina/R2.trim.fastq.gz | \
  samtools sort -@ 16 -o assembly/illumina.sorted.bam
samtools index assembly/illumina.sorted.bam

pilon --genome assembly/flye/assembly.fasta \
      --frags assembly/illumina.sorted.bam \
      --outdir assembly/pilon_r1 --threads 16

quast -o assembly/quast assembly/pilon_r1/pilon.fasta
busco -i assembly/pilon_r1/pilon.fasta -l fungi_odb10 -m genome -o assembly/busco

```
Expected: few contigs, high N50, BUSCO completeness >90%.
Theory link: single-copy orthologs benchmark assembly completeness.
Galaxy alternative: Flye, Pilon, QUAST, BUSCO wrappers.

## 5) Annotation

Concept: MAKER3 integrates ab initio predictors (AUGUSTUS, GeneMark-ES), RNA-Seq, and protein homology. dbCAN3 annotates CAZymes.

Procedure (mock):

```bash
maker -CTL
# edit maker_opts.ctl for genome, proteins, transcripts
maker
gff3_merge -d *.maker.output/*/*_master_datastore_index.log > maker_merged.gff
fasta_merge -d *.maker.output/*/*_master_datastore_index.log > maker_merged.fasta

run_dbcan maker_merged.fasta protein --out_dir dbcan_results --cpu 16
```

Expected: 10–15k gene models, enriched in esterases, cutinases, lipases.
Theory link: integration reduces false positives.
Galaxy alternative: AUGUSTUS, InterProScan, dbCAN.

## 6) RNA-Seq Quantification & DEG Analysis

Concept: Salmon quantifies reads; DESeq2 identifies DEGs.

Procedure (mock):

```bash
salmon index -t annotation/maker/maker_merged.fasta -i rna/index --type quasi -k 31

salmon quant -i rna/index -l A \
  -1 data/rna/sample_R1.fastq.gz -2 data/rna/sample_R2.fastq.gz \
  -p 16 -o rna/quants/sample

# R/DESeq2
R --vanilla <<'RSCRIPT'
library(tximport); library(DESeq2)
samples <- read.csv("metadata/rna_samples.csv")
files <- file.path("rna/quants", samples$sample, "quant.sf")
names(files) <- samples$sample
txi <- tximport(files, type="salmon", txOut=TRUE)
dds <- DESeqDataSetFromTximport(txi, colData=samples, design=~condition)
dds <- DESeq(dds)
res <- lfcShrink(dds, coef=2, type="apeglm")
write.csv(as.data.frame(res), file="rna/deseq2_results.csv")
RSCRIPT
```

Expected: ~300 DEGs; hydrolases strongly induced.
Theory link: DEGs link genomic content to functional activity.
Galaxy alternative: RNA-Seq → DESeq2 workflows.

## 7) Comparative Genomics

Concept: OrthoFinder identifies orthogroups; MCScanX analyzes synteny; Mauve detects rearrangements.

Procedure (mock):

```bash
orthofinder -f comparative/genomes/ -t 32 -a 32 -o comparative/orthofinder

makeblastdb -in proteomeA.faa -dbtype prot
blastp -query proteomeB.faa -db proteomeA.faa -outfmt 6 -evalue 1e-5 -num_threads 32 > A_vs_B.blast
MCScanX A_B
```

Expected: unique orthogroups in degraders; conserved synteny with local expansions.
Theory link: orthology and synteny reveal adaptive signatures.
Galaxy alternative: OrthoFinder, synteny modules.

## 8) Phylogenomics

Concept: Build phylogenies from single-copy orthologs using ML and Bayesian approaches.

Procedure (mock):

```bash
mafft --thread 32 input_single_copy.faa > aligned.faa
Gblocks aligned.faa -t=p -b5=h
raxmlHPC -T 32 -s concatenated.phy -n tree -m PROTGAMMAJTT -p 12345 -# 100

# MrBayes
mb <<'MB'
execute tree.nex
lset rates=gamma aamodelpr=fixed(jones)
mcmc ngen=2000000 samplefreq=1000 nchains=4
sump; sumt;
MB
```

Expected: degraders cluster into a distinct clade.
Theory link: single-copy orthologs provide robust phylogenetic signals.
Galaxy alternative: MAFFT, RAxML, MrBayes.

## 9) Integration

Expected Narrative:

Genome assemblies → provide the toolbox

Annotation → reveals candidate enzymes

Transcriptomics → shows inducible expression

Comparative genomics → highlights lineage-specific traits

Phylogenomics → places isolates in evolutionary context

Theory link: integrated pipelines demonstrate how genomics connects molecular data to applied biotechnology.
Galaxy alternative: chain workflows end-to-end.

## 10) HPC Templates

Concept: Submit heavy jobs on SLURM clusters.

```bash
#SBATCH -J flye
#SBATCH -c 32
#SBATCH --mem=180G
#SBATCH -t 48:00:00
module load mamba && mamba activate agproj
flye --nano-raw data/nanopore/reads.fastq --out-dir assembly/flye --threads 32

#SBATCH -J orthofinder
#SBATCH -c 48
#SBATCH --mem=256G
#SBATCH -t 72:00:00
module load mamba && mamba activate agproj
orthofinder -f comparative/genomes/ -t 48 -a 48 -o comparative/orthofinder
```

## 11) Galaxy Shortcuts

QC: FastQC → MultiQC

Trimming: Fastp, Cutadapt

Assembly: Flye + Pilon

QC: QUAST, BUSCO

Annotation: AUGUSTUS, InterProScan, dbCAN

RNA-Seq: Salmon → DESeq2

Comparative: OrthoFinder, MAFFT, RAxML, synteny

## 12) Final Notes

All commands are mock examples; adapt to your cluster, dataset, and resources.

BUSCO requires lineage DBs (e.g., fungi_odb10).

GeneMark-ES requires a license.

This pipeline mirrors Applied Genomics course modules: sampling, QC, hybrid assembly, annotation, transcriptomics, comparative genomics, phylogenomics, and integration.
