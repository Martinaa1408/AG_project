# Applied Genomics – Slide Recap and Concept

Full slide deck: [Google Drive](https://drive.google.com/file/d/1tV58Ldxbase2jpN9sgMsIyVQvUfpziG5/view?usp=sharing)

---

## 📑 Repository Structure

- **mock exam/**  
  Contains mock exams and practice questions.  
  --> Useful to test your understanding and prepare for the written exam.

- **notebook/**  
  Collection of course notebooks
  --> Base material to review concepts.

- **project/**  
  Final project folder  
  --> Includes presentation and some details.

---

**Genetics** is the discipline that studies **biological diversity at different levels: genomes, populations, individuals and organisms**.
**Diversity** is the **driving-force** of genetics.

---

## 1. Branches of Genetics

## **Classical genetics (transmission/formal)**: Mendelian inheritance, pedigree analysis, laws of segregation & independent assortment.

  -->heredity and segregation of traits.
  
     - focuses on the **principles of heredity** and how traits are transmitted.
     - **chromosomes** are the physical carriers of heredity, containing multiple genes.
     - **crossing over** during meiosis reshuffles paired chromosomes; recombination frequency depends on distance betweeen loci.
     - **first genetic maps** were developed in Drosophila later replaced by genome sequencing.
     - loci were studied by phenotype changes-->Only **loci with different alleles** could be detected in classical genetics.
       
     - Mendel formulated the **laws of heredity** from Pisum sativum experiments-->
              -**Genotype** (allele combination) vs **Phenotype** (observable traits: genotype + environment)
              -**Dominant** allele masks recessive; **Codominance**-->intermediate phenotypes possible
              -**1) Law of segregation**: each gamete carries one allele form each locus; classic 3:1 phenotype ratio in F2 generation.
              -**2) Law of Independent Assortment**: alleles for different traits assort independently if loci are on different chromosomes; 9:3:3:1 in dihybrid crosses.
       
     - **Pedigrees**: graphical representation of inheritance across generations.
         -Symbols: o female and square male rombous unknown sex; empty=healthy filled=affected and half-filled=carrier
         -**PLINK** is the software for genotype-phenotype data management;
            input: text file (one raw per sample and multiple fields); columns=family, individual ID, parents, sex, phenptype.
            facilitates high-throughput genotyping analysis and genome-wide association studies (GWAS).
     - **cytogenetics**: visulization of chromosomes during metaphase (condensed); involves staining techniques to identify structure and number.
     - **chromosomal map**: represent gene positions and distances based on recombination frequency; closely linked genes-->low recombination.
     - **Sex chromosomes** and determination: in mammals XX is the female and XY the male, PARR (pseudoautosomal regions) allow delimited recombination between X and Y.
         in other system X0 (insect): XX female and X male; ZW (birds): ZZ male and ZW female; haplo-dyploid (bees): fertilized diploid female, unfertilized haploid male;           temperature-dependent sex in reptiles.

<img width="374" height="364" alt="Screenshot 2025-08-29 153600" src="https://github.com/user-attachments/assets/e3338577-e990-4e36-afa0-806ab9d49a26" />

       
---

## **Molecular genetics**: DNA as genetic material, PCR, Sanger sequencing.
  
  -->focuses on the chemical nature of genes, their structure, function, and how genetic information is encoded, replicated, expressed and regulated.
  
     - **DNA replication, transcription, and translation**.
     - **gene regulation**.
     - **post-transcriptional and post-translational processes**.
     - the technical foundations of molecular genetics include: **recombinant DNA technology**, **DNA sequencing**, **library construction**, **PCR amplification**,               ***hybridization techniques**, **gel electrophoresis**.
     - **SANGER SEQUENCING: first generation sequencing, an early but fundamental sequencing method**.
           Principle:
            - based on controlled DNA synthesis with chain-terminating nucleotides.
            - requires: DNA template (to be sequenced), a primer (defines start point); dNTPs normal nucleotides, ddNTPs (modified nucleotides, terminate synthesis) and                            DNA polymerase enzyme
           Process:
           1- run four reactions (one per nucleotide);
           2- incorporate normal or modified nucleotides randomly during synthesis;
           3- modified nucletoides block extension, generating DNA fragments of different lengths.
           4- fragments are separated by electrophoresis which sorts by size;
           5- detection is based on light or color emitted by the final ddNTP;
           6- by analyzing fragment lengths and terminal nucleotides the DNA sequence is reconstructed.
       -->key features: uses linear amplification millions of fragments of varying lengths, highly accurate but low throughput compared to modern methods, provided the               foundation for today's sequencing technologies.
   <img width="376" height="308" alt="Screenshot 2025-08-29 153827" src="https://github.com/user-attachments/assets/7dfb8363-8369-406a-8f9c-abada9e001d9" />
    
     - **PCR (polymerase chain amplification)**: method for amplification. Its goal is to take a specific tiny segment of DNA and make billions copies of it. This process           is exponential amplification, each cycle double the number of copies.
           Process:
             - denaturation --> annealing --> extension
           Principle:
             - requires: template DNA, primers, Taq DNA polymerase, dNTPs and buffer
       
     --> link: PCR is first used to isolate and massively amplify the specific DNA region of interest, then Sanger is used to determine the precise genetic code of that             amplified PCR product.
<img width="350" height="300" alt="Screenshot 2025-08-29 153819" src="https://github.com/user-attachments/assets/16f43640-5a2a-4f91-ab4f-1479ec366ede" />

---

## **Population genetics**: Hardy–Weinberg equilibrium, allele/genotype frequencies, LD, inbreeding, ROH.
  
  -->variability within and across populations, infact studies the genetic composition of populations (group of individual of the same species) and howthis composition          changes over space and time

     - focus on **allele frequencies** in populations.
     - allele frequencies can be used to infer **genotype frequencies** under certain assumptions.
     - **HWE (hardy weinberg equilibrium)**: a central model (null model) in population genetics to predict expected genotype frequencies from allele frequencies.
       Assumptions: the population must satisfy the **following conditions**:
         - diploid organisms
         - only sexual reproduction occurs
         - non overlapping generations
         - random mating
         - infinite population size
         - equal allele frequencies in sexes
         - no migration, selection, mutation or gene flow
      If all assumptions hold, allele and genotype frequencies remain constant across generations.
       **Allele frequencies = Number of copies of a specific allele/total number of alleles at that locus in the population**
       **Genotype frequencies = Number of individuals with a given genotype/Total number of individuals in the population**
      Ex. biallelic locus: let p= frequency of the dominant allele and let q= frequency of the recessive allele, since only two alleles exist: p + q = 1.
      The expected genotype frequencies follow: **(p+q)^2= p^2 + 2pq + q^2 = 1**.
        where p^2 = homozygous dominant (AA), 2pq = Aa heterozygous (Aa), q^2 = homozygous recessive (aa).
      Interpretation: if the observed population frequencies differ from HWE expectations,this indicates that evolutionary forces may be acting; the model is not a               description of real populations but a reference baseline to study deviations.

     - **LD (linkage disequilibrium)**: non-random association of alleles at different loci in a population. Stronger LD when loci are physically closer tend to be                inherited together.
       Influenced by:
        - recombination rate (low recombination=high LD).
        - population size (small effective size= high LD).
        - inbreeding and bottlenceks (increase LD).
     - **Inbreeding coefficients**: Fped (pedigree-based probability that two alleles are identical by descent); FROH (fraction of genome covered by ROH); Fis (inbreeding         coefficient measuring excess homozygosity relative to HWE.
     - LD (r², D′): non-random association of alleles at loci.
     - **ROH (runs of homozigosity)**: longs continuous stretches of homozygous genotypes in the genome; indicate inbreeding both chromosomes inherited from a common              ancestor. Long ROH = recent inbreeding and short ROH = ancient inbreeding.
     - **centiMorgan (cM)**: genetic distance unit based on recombination frequency. 1 cM = 1% recombination (chance of crossover per generation). Loci close together show        high LD.
     - **F_ST**= measures genetic differentiation between populations
       Var(p)/p(1-p) ranges: 0 no differentiation and 1 populations are fixed for different alleles.

LINK-->HWE is the theoretical baseline for LD, ROH, and inbreeding measures; LD is shaped by recombination measured in cM and mantained by inbreeding or drift; cM quantifies recombination which breaks down LD over time; ROH provide a genomic measure of inbreeding; all capture inbreeding but from different sources (pedigree, genomic data, frequency-based); F_ST extends the concept of inbreeding and HWE deviations to the population level.

[HWE] baseline equilibrium-->(deviations)-->[inbreeding]--[ROH] genomic signal-->(quantified by Fped,FROH,FIS)-->(reduces effective recombination)-->[LD high]-->(linked to)-->[cM]-->(across populations)-->[F_st] differentiation between populations

Population applications:  
- **Genotyping (SNP arrays, GBS, RAD)** → captures genome-wide variants; essential for diversity studies, structure, and as input for GWAS. Arrays = known variants (cheap, reproducible), GBS/RAD = discover SNPs in non-model species.  
- **Structural variants (CNVs, inversions, translocations)** → large-scale changes often affect gene dosage or regulation; critical in cancer genomics, evolution, and adaptation. Detected with arrays, short reads (depth, discordant pairs), or more effectively with long reads.  
- **GWAS** → links genotypes to phenotypes without prior hypotheses. Pipeline: QC → correction for structure (PCA, kinship) → association → multiple testing correction → visualization (Manhattan). Strength: unbiased trait discovery; Limits: rare variants hard to detect, risk of stratification bias.  
- **MAS / MAB** → bridges population data to breeding. MAS selects individuals with markers linked to traits; MAB/genomic selection uses genome-wide SNPs + predictive models to estimate breeding values (GEBVs), speeding up improvement of polygenic traits.  
- **Conceptual link to population genomics** → LD underlies SNP associations; ROH and inbreeding shape background signals; F_ST reflects differentiation and thus resolution of associations across populations.  

---
       
## **Quantitative genetics**: phenotype = genetics + environment, variance decomposition, heritability, QTL.
  
-->polygenic traits, studies traits influenced by many genes and environmental factors (ex. stature, eye color, running ability, survival to adulthood)
     
    Traits can be continuous (height), discrete (litter size), or binary (survival yes/no).
    -key concept: **phenotype (P) = Genotype effect (G) + Environment effect (E)**
      component: **genotype effect** = **addictive effect (A)** sum of contributions of alleles across loci; **dominance effect (D)** interaction between alleles at the                     same locus; **interaction effect or epistatic (I)** gene-gene interactions.
                 **environment effect** = permanent environmental effect and temporary environmental effect.
     ---> thus: Var (P) = Var (G) + Var (E)
    - **Heritability** is the proportion of phenotypic variance explained by genetic variance. (h^2=Var(G)/Var(P))
      low heritability: < 0.1
      medium heritability: 0.1 - 0.4
      high geritability: > 0.4
      To estimate heritability, genetic realtionships between individuals must be considered.
    - Application: quantitative genetics is essential to understand variation and covariation among relatives in natural and managed populations, study the dynamics of           evolutionary change, improve animals and plants and investigate complex disease in humans.
      
---
     
## **Genomics**: 

whole-genome study, WGS, metadata, multi-omics.
-->what is the **genome**? the entire genetic content of an organism (DNA sequence).
-->what is **genomics**? the study of all genes in an organism and their interactions. Is a multidisciplinary field combining: genetics, molecular biology, robotics and computing.

**1990-2001 Human Genome Project (HGP)**: sequenced 3 billion base pairs, high quality (<1 error per 10k bases)

**Omics expansion**: functional genomics (gene expression and function), transcriptomics (RNA profiles), proteomics (protein complement), metabolomics (metabolites), phenomics (phenotype data)

**diversity of genomes**: 3 domains of life: **bacteria** (simple but share some genes with humans), **eukaryota** (complex genomes-->nuclear + organellar), **archea** (extreme environments, thermophiles, helophiles and methogens).

**types of genomes**: viruses (DNA or RNA genomes), prokaryotes (circular genomes + plasmids no reference genome concept), eukaryotes (nuclear genome the main DNA, mitochondrial 60 kb, chloroplast).

**data and databases**: genebank repository of dna sequences with metadata (species, source, annotations); WGS project randomly sequence dna fragments and assemble later; comparative genomics compare annotated and unannotated genomes to infer gene functions.

**big data in genomics**: NGS has made genomic data astronomical: Terabyte (TB) = 10^12 bytes, Petabytes (PB) = 10^15 bytes, Exabytes (EB) = 10^18 bytes, Zettabytes (ZB) = 10^21 bytes (1 million TB)
Challenges: Storage, Transfer, Analyses.

**metadata**=contextual information accompanying genomic data.

---

## Branches of Genetics – Conceptual Links

| Branch                  | Core Idea / Remember Link                                                                                                      |
|--------------------------|-------------------------------------------------------------------------------------------------------------------------------|
| **Classical genetics**   | **Mendel → Laws (Segregation, Independent Assortment)** → inheritance of traits → chromosomes as carriers → crossing-over → genetic maps → pedigrees → sex determination systems |
| **Molecular genetics**   | **DNA → Information flow (Replication → Transcription → Translation)** → gene regulation → PCR amplifies → Sanger reads sequence → foundation of modern sequencing & recombinant DNA |
| **Population genetics**  | **HWE (p² + 2pq + q² = 1) → baseline equilibrium** → deviations = evolutionary forces → LD (non-random alleles) → recombination (cM) → inbreeding (Fped, FROH, Fis) → ROH = genomic inbreeding signal → F_ST = differentiation among populations |
| **Quantitative genetics**| **Phenotype (P) = Genotype (G) + Environment (E)** → variance decomposition (A+D+I+E) → heritability (h² = VarG/VarP) → QTL detect loci → polygenic traits (continuous, discrete, binary) → applied to breeding & complex disease |
| **Genomics**             | **Genome = all DNA** → Human Genome Project (HGP, 1990–2001) → sequencing technologies (NGS) → big data challenges (TB → ZB) → multi-omics (functional, transcriptomics, proteomics, metabolomics, phenomics) → comparative genomics & databases (GenBank, WGS) |

---

## 2. Sequencing Technologies

## **NGS next generation sequencing**
is the main technology in modern genomics, enabling the production and analysis of massive amounts of sequencing data; requires consideration of error rates throughput costs and timelines and provides unbiased sequencing of millions of DNA fragments simultaneously.
  -many DNA targets sequenced in parallel
  -millions of fragments analyzed in a single run
  -uses clonal PCR amplification to boost signal detection
  -sequences are unknown in advance identified only after computational analysis.
  
## **Sanger (1st generation)**
  -accurate, low throughput; one target DNA-->one PCR-->one sequencing reaction.
  -each capillary= one lane max 58 sequences/run
  -targeted: sequence is know in advance (via specific primers)
  -expensive, accurate, low throughput.
  
<img width="316" height="395" alt="image" src="https://github.com/user-attachments/assets/7cba553b-ca39-454d-b0e2-ffe2741f5832" />

  
LINK: from 'we know what we sequence (SANGER)' to 'we discover what we sequenced (NGS)'.

**data production**: modern sequencers up to 48 samples/run; run duration is 1.5 hours and read length 1.5 kb per sample.
Every sequencing experiments involves 3 stages: Planning, data production and analyses.
In sanger era the data production was the most expensive and limiting factor and few sequences produced long timelines and limited insights
In present NGS data production is cheap and fast planning is crucial to optimize study deisgn and analysis is the bottleneck.

Ex. HGP 100,000,000 per genome today--> <1000 per genome, cost per 1 Mb DNA <0,01 
-->NGS has dramatically outspaced MOOre's law in reducing costs.

<img width="355" height="360" alt="image" src="https://github.com/user-attachments/assets/cfa7e9fe-4b68-4071-bc29-bfe6b73ad0eb" />


| Aspect                  | **Sanger Sequencing (1st Gen)**                              | **NGS (Next Generation Sequencing)**                            |
|--------------------------|-------------------------------------------------------------|-----------------------------------------------------------------|
| **Concept**             | *We know what we sequence* (target defined by primers)       | *We discover what we sequenced* (massively parallel, unbiased)  |
| **Throughput**          | Low: one DNA → one PCR → one sequencing reaction; max ~96 capillaries/run (~58 seqs/run) | High: millions of fragments sequenced in parallel in a single run |
| **Input / Target**      | Specific known sequence (targeted, primer-based)             | Random, genome-wide, no prior knowledge needed                  |
| **Signal amplification**| PCR product sequenced directly                               | Clonal amplification (emulsion PCR, bridge PCR) → boosts signal |
| **Data production**     | Few sequences, expensive, time-consuming                     | Dozens of samples/run, fast (~1.5h), cost-efficient, high yield |
| **Read length**         | ~800–1000 bp per read                                        | Short reads (100–400 bp) or long reads (PacBio/Nanopore)        |
| **Error rate**          | Very low (high accuracy)                                     | Higher error rate (depends on platform), computational correction needed |
| **Main limitation**     | Costly, low throughput, limited discovery power              | Analysis bottleneck: huge data → requires planning & computation |
| **Era bottleneck**      | Sanger: **data production** (slow, costly, few insights)     | NGS: **analysis** (planning crucial, data cheap & fast)         |

---

## **Ion Torrent sequencing system**: detects H+ release, ionograms.
  Ion torrent (now Thermo Fisher Ion S5 system) is a **short-read NGS technology based on semiconductor sequencing**.
  Instead of using fluorescence or optics, the system detects pH changes caused by nucleotide incorporation.
  Sequencing occurs in tiny wells on a silicon ship, each acting as an independent reaction chamber.
  Each well= one DNA fragment per thousands of clonal copies (to amplify signal)

  -->**chips**: the chip is the core of the ion torrent system; each well = one sequencing reaction with a pH sensor to detect ions; dna fragments are loaded into wells         and sequenced in parallel; different chip models= different throughput scalability (510 2,5 Gb output 2.3 million reads).
  -->**sequencing principle**: dna fragments are engineered with universal primer sites, one nucleotide type is flowed across the chip at a time. If complementary the           polymerase adds it to the DNA strand releasing H+ ions detected as pH change and pyrophosphate, each pH change=incorporation event converted into digital signal.           Signal intensity is proportional to the number of identical nucleotides added (homopolymers)
  -->**workflow**: DNA fragmentation (sonication or enzymatic)-->library prep (adapters + barcodes)-->emulsion PCR(clonal amplification)-->chip loading (dna loaded beads        into wells)-->sequencing(ion detection via pH sensors)-->data output (torrent suite software for qc and read statistics).
  -->**library and template prep**: dna fragmentation require sonication (physical) or enzymatic digestion, then library preparation require fragments engineered with           adapters + universal primer sites, barcoding allows multiplexing samples from multiple individuals and clonal amplification performed by emulsion PCR (one fragment         DNA per micro-droplet), creates thousand of copies of the same DNA fragment and the goal is 1 dna fragment per sphere avoids polyclonals (mixed signals).
  -->**flow and sequencing**: sequencing occus by programmed flows of nucleotides; homopolymer regions (AAA, TTT) signal proportionality problem higher error rate; output       is displayed in **ionograms** similar to elctropherograms.
  -->**strengths**: fast and relatively cheap (no fluorescently modified nucleotides and block), simple chemistry, highly scalable with different chips and support             barcoding (cost-efficient multiplexing)
  -->**limitation**: error-prone in homopolymers regions (cannot distinguish 3 vs 4 identical bases), require careful DNA quantification to avoid polyclonals (multiple         template per sphere and empty wells, chips are single-use, read-length limited to 200-400 bp.

<img width="394" height="301" alt="Screenshot 2025-08-29 154419" src="https://github.com/user-attachments/assets/732f4722-5325-4a6a-b51c-95298512f151" />

| Aspect                  | Key Points / Keywords                                                                 |
|--------------------------|---------------------------------------------------------------------------------------|
| **Technology**           | Short-read NGS, **Sequencing by Synthesis (SBS)** → detects **H⁺ ions (pH change)** instead of fluorescence/optics |
| **Chip**                 | Silicon chip with wells = independent reaction chambers; each well = 1 DNA fragment (thousands of clonal copies); pH sensor detects incorporation; different chip models = scalable throughput (e.g. 510 chip ≈ 2.5 Gb, ~2.3M reads) |
| **Sequencing principle** | DNA fragments with universal primers; nucleotides flowed sequentially; if incorporated → polymerase adds base → H⁺ released → pH change detected; **signal intensity ∝ # of incorporated bases (homopolymers issue)**; output = **ionograms** |
| **Workflow**             | DNA fragmentation (sonication/enzymatic) → library prep (adapters + barcodes) → **emulsion PCR (clonal amplification)** → chip loading (DNA beads into wells) → sequencing (pH detection, SBS method) → data analysis (Torrent Suite software) |
| **Library/Template prep**| DNA fragments + adapters + universal primer sites; barcoding for multiplexing; clonal amplification by emulsion PCR (1 fragment per droplet → 1 bead with thousands of identical copies); avoid polyclonals/empty wells |
| **Flow & sequencing**    | Programmed flows of single nucleotide types → incorporation → H⁺ release → signal; **homopolymer regions (AAA/TTT)** cause signal proportionality error; output visualized as ionograms (similar to electropherograms) |
| **Strengths**            | Fast, relatively cheap, simple chemistry (no fluorescent dyes/blockers), scalable with different chips, barcoding enables multiplexing |
| **Limitations**          | Homopolymer error (3 vs 4 bases indistinguishable), requires accurate DNA quantification, risk of polyclonals, single-use chips, short read length (200–400 bp) |

-->**Electropherogram**= fluorescence colors → high accuracy → low throughput
-->**Ionogram** = pH intensity → high speed/parallelism → homopolymer error risk

---

## **Illumina SBS**: short, accurate reads, cluster generation, paired-end.
  Illumina is the **leading NGS platform** in genomics today. It is based on **Sequencing by Synthesis (SBS)** where nucleotides are incorporated one at a time detected by   fluorescence and then reset for the next cycle.
  -produces short-reads (100-300 bp)
  -extremely high accuracy (Q30-99.9%)
  -scalable with different machines and flow cells
  -->**workflow**: 1-library prep: DNA fragmentation by physical sonification or enzymatic digestion; adapters are added on DNA fragments enable binding to flow cell            oligos, contain barcodes for multiplexing and provide primer binding sites.
     2-cluster generation (bridge amplification): performed on a flow cell (glass slide coated with oligonucleotides) each fragment binds the surface and forms a bridge,        PCR amplification occurs creating thousand of clonal copiesd of each fragments and each cluster = one DNA fragment amplified, ordered clusters in modern flow cells --      >higher accuracy and resolution.
     3-sequencing by synthesis: reversible terminator nucleotides each base has a fluorescent label + chemical blocking group, only one nucleotide is incorporated per           cycle. Workflow per cycle: -added 4 fluorescently labeled nucleotides; -DNA polymerase incorporates the complementary base; imaging system captures the fluorescent         color of each cluster; -blocking group is removed and synthesis can continue; -cycle repeats. Output: fluorescence signals translated into read sequence.
  -->**chemistry and improvements**: 4 channel SBS classic method-->4 fluorochromes one per nucleotide, 2 channel SBS uses only 2 dyes for 4 bases reduced cost/time, 1          channel SBS uses one dye with chemical tricks (on/off states) to distinguish nucleotides. The advantages of SBS chemistry is to prevent homopolymers errors, high           precision because reaction is stopped and read at each cyle. Trade-off: chemistry is more expansive due to modified nucleotides, sequencing is slower because of stop-      read-reset-cycles.
  -->**accuracy and read length**: error rate extremely low, gold standard for WSG, exome sequencing, RNA-Seq and small variant detection (SNPs, indel).
  -->**advantages**: very high accuracy, wide range of applications, scalable throughput (different flow cells and instruments), multiplexing via barcodes.
  -->**limitationss**: more expansive chemistry (modified nucleotides), short-reads compared to long-read technologies, slower per cycle (due to imaging and chemical reset      steps).
  -->fluorescent label attach to the base, and blocking group in the 3'.

<img width="332" height="322" alt="Screenshot 2025-08-29 154648" src="https://github.com/user-attachments/assets/b764f4e1-79a5-4b04-b97a-161f69245ce6" />

<img width="374" height="323" alt="Screenshot 2025-08-29 154732" src="https://github.com/user-attachments/assets/d40825e2-6948-491c-b1cb-43a1fa85806e" />



| Aspect                  | Key Points / Keywords                                                                 |
|--------------------------|---------------------------------------------------------------------------------------|
| **Technology**           | Leading NGS platform, **Sequencing by Synthesis (SBS)** → nucleotides incorporated one at a time, detected by **fluorescence**, reset each cycle |
| **Read type**            | Short reads (100–300 bp), **paired-end** option, Q30 ≥ 99.9% (extremely accurate)     |
| **Workflow**             | 1. **Library prep** → DNA fragmentation (sonication/enzymatic) + adapters (barcodes, primer sites) <br> 2. **Cluster generation** → bridge amplification on flow cell, clonal copies per cluster <br> 3. **SBS sequencing** → reversible terminator nucleotides, fluorescence detection per cycle, blocking group removed, repeat |
| **Chemistry**            | Classic: **4-channel SBS** (4 dyes, one/base) <br> 2-channel SBS: 2 dyes → 4 bases <br> 1-channel SBS: 1 dye with on/off chemical states <br> **Advantage**: prevents homopolymer errors, precise base-by-base reading <br> **Trade-off**: chemistry expensive, slower cycles (stop → read → reset) |
| **Accuracy & length**    | Very low error rate, gold standard for WGS, exome, RNA-Seq, SNP/indel detection       |
| **Advantages**           | Ultra-high accuracy, wide applications, scalable throughput (machines & flow cells), supports multiplexing via barcodes |
| **Limitations**          | Expensive chemistry (modified nucleotides), short reads (assembly harder than long-read tech), slower cycles due to imaging & reset |

-->**Illumina Paired End Sequencing**: Illumina sequencing is usually short-read sequencing. To overcome the limitation of sequencing only small fragments paired-end sequencing allow reading both ends of a DNA fragment.
Ex. 1000 bp fragmented DNA illumina paired end sequence first 200 bp and last 200 bp the middle 600 bp remain unsequenced directly but are physically linked.
2 short reads separated by an insert size (unsequenced gap).
Applications: genome assembly (de novo, resequencing), structural variation (reveal translocations, inversion, gene fusions), transcriptomics (splicing events),improved reconstruction of missing sequences.
-->**Mate-Pair Sequencing (Illumina special library prep)** → Designed for longer insert sizes (2–20 kb). DNA is circularized, fragmented, and sequenced at the junctions. As a result, the two reads come from far apart in the original genome (kb distance), not just hundreds of bp like paired-end.
Applications: de novo genome assembly scaffolding, resolving repetitive regions, detecting large structural rearrangements, mapping across long gaps.
-->Link between them:
Paired-end = short insert (~200–600 bp) → precise, high-resolution, local context.
Mate-pair = long insert (2–20 kb) → long-range linking information, helps scaffold assemblies and detect large structural variation.

<img width="321" height="398" alt="Screenshot 2025-08-29 155526" src="https://github.com/user-attachments/assets/dd694247-e4e6-4baf-996f-6727bd7f9d46" />

---

## **454 Roche & ABI SOLiD**: early NGS, discontinued.
  **ABI SOLiD**: sequencing by ligation system, It identifies the sequence by matching and ligating pre-made probes with ligase.
    Its key defining feature was its use of DNA ligase instead od DNA polymerase to determine the sequence, and its unique two-base encoding system for color space. The        throughput is impressive but have complex data analysis and short-reads (50-75 bp).
    -->Step: emPCR; bind beads to a glass slide; hybridize and ligate fluorescent probes; detect color from ligation.
  **454 Roche**: SBS method, It build the DNA strand by adding nucleotides with polymerase.
    Ion Torrent is to faster and significantly cheaper as it eliminated the need for expansive enzymatic reagents and complex optical imaging systems.
    -->Step: emPCR; load beads into Picotiterplate wells; flow nucleotides; detect light from incorporation.
    -->the detection method is pyrosequencing PPi (a biochemical reaction that produces light)


| Platform        | Principle / Chemistry                                                                 | Workflow Steps                                                                                          | Features & Read Length           | Pros                                                                 | Cons                                                                 |
|-----------------|---------------------------------------------------------------------------------------|---------------------------------------------------------------------------------------------------------|----------------------------------|----------------------------------------------------------------------|----------------------------------------------------------------------|
| **ABI SOLiD**   | **Sequencing by Ligation**: uses **DNA ligase** + fluorescent probes; unique **two-base encoding (color-space)** | 1. emPCR (clonal amplification on beads) <br> 2. Beads immobilized on slide <br> 3. Fluorescent probes hybridized & ligated <br> 4. Color signals detected | Short reads **50–75 bp**         | High throughput, ligation-based accuracy (error-correction via color-space) | Complex data analysis (color-space), short reads, platform discontinued |
| **454 Roche**   | **Sequencing by Synthesis (SBS)** via **pyrosequencing** (detects light from PPi release) | 1. emPCR (DNA on beads) <br> 2. Load beads into **PicoTiterPlate wells** <br> 3. Flow nucleotides sequentially <br> 4. Incorporation → PPi release → luciferase reaction → light detected | Reads up to **400–700 bp**       | Longer reads than Illumina at the time; fast compared to Sanger; scalable with picotiterplate | Expensive reagents (luciferase/enzymes), homopolymer errors, high cost per base, discontinued |

<img width="324" height="374" alt="Screenshot 2025-08-29 154640" src="https://github.com/user-attachments/assets/af9b64f4-8fc5-4b3f-b110-daed4603d242" />
<img width="350" height="397" alt="Screenshot 2025-08-29 154615" src="https://github.com/user-attachments/assets/012f65fc-3cac-4d31-be22-ecbe900f5662" />

---

**Long-read sequencing**: while Illumina dominates short-read sequencing two main technologies provide long reads (up to tens of kb or more): ONT and PAcBio.
Long read sequencing helps resolve: repetitive regions, structural variations, complex genome assemblies, phasing of haplotypes.

## **PacBio (SMRT, HiFi)**: long accurate reads.
  -->Principle: a single DNA molecule is immobilized in a well (Zero-Mode Waveguide, ZMW). DNA polymerase incorporates fluorescently labeled nucleotides. Each                   incorporation event emits light captured by a camera. Interpulse duration, time between signals, helps determine the base.
  -->**CCS Circular Consensus Sequencing**: DNA fragments are circularized SMRTbell templates; sequenced multiple times as polymerase loops around the circle; errors are        random corrected by overlapping multiple passes.
  -->Features: read length typically 10-25 kb, can reach 50 kb. Error correction: high accuracy with CCS 8HiFi reads). Bias reduction: no issues with GC-rich regions.           Phasing: can separate maternal and paternal haplotypes.
  -->Pros: accurate long reads (HIFi reads with Q30+); excellent for strctural variants and complex genomes; NO GC-bias; phasing of alleles possible.
  -->Cons: high cost (3-4x Illumina per genome); lower throughput than short-read systems; library prep is complex and requires high-quality non degraded DNA; larger lab-       based instruments.
  -->fluorescent label on the phosphate chain (not on the base)
  
 <img width="310" height="393" alt="Screenshot 2025-08-29 154902" src="https://github.com/user-attachments/assets/7bd63742-f894-49e2-953b-743424791ec2" />
<img width="362" height="348" alt="Screenshot 2025-08-29 154931" src="https://github.com/user-attachments/assets/508c3567-401f-4c0f-b500-a187a8cf2341" />


| Aspect                  | Key Points / Keywords                                                                 |
|--------------------------|---------------------------------------------------------------------------------------|
| **Principle**            | **Single-Molecule Real-Time (SMRT)**: DNA immobilized in **Zero-Mode Waveguide (ZMW)** wells; polymerase incorporates fluorescent nucleotides → light pulses detected; interpulse duration aids base calling |
| **CCS / HiFi**           | **Circular Consensus Sequencing (CCS)**: DNA circularized into SMRTbell templates → polymerase makes multiple passes → random errors corrected → **HiFi reads (Q30+, high accuracy)** |
| **Read length**          | Typical: **10–25 kb**; up to **50 kb** possible                                      |
| **Unique features**      | Random error correction with CCS; no GC-bias; **haplotype phasing** (maternal vs paternal separation) |
| **Pros**                 | Accurate **long reads** (HiFi Q30+); excellent for structural variants & complex genomes; unbiased (GC-rich ok); haplotype phasing |
| **Cons**                 | High cost (3–4× Illumina); lower throughput; requires high-quality DNA; complex library prep; large lab instruments |



## **Nanopore (MinION, PromethION)**: ultra-long reads, real-time.
  -->Principle: sequencing native DNA/RNA strands (no synthesis step). A single stranded DNA passes through a nanopore embedded in a membrane. Electric current runs across      the pore each nucleotide or combination of bases causes a specific voltage disruption. Volatge shifts are translated into base calls.
  -->Features: read length up to 10-20 kb in practice, theoretical maximum >1 Mb. Portability devices like MinION are samll, USB-powered and usable in the field. Direct         sequencing of RNA (no cDNA needed).
  -->Pros: very cheap instrument, portable usable on-site, ultra long reads possible, sequencing of RNA directly
  -->cons: high error rate (5% historically 20%), throughput less stable, library prep still required and not suitable for degraded samples.
  
  <img width="387" height="342" alt="Screenshot 2025-08-29 154839" src="https://github.com/user-attachments/assets/d0b1a290-1deb-4b60-af0c-e16bde077bf7" />
<img width="318" height="381" alt="Screenshot 2025-08-29 154823" src="https://github.com/user-attachments/assets/b2137387-08b0-4232-8eee-44d6bfdee37d" />


| Aspect                  | Key Points / Keywords                                                                 |
|--------------------------|---------------------------------------------------------------------------------------|
| **Principle**            | Native **DNA/RNA strands** pass through nanopore in a membrane; electric current across pore → each nucleotide causes **voltage disruption** → signal translated into bases |
| **Read length**          | Practical: **10–20 kb**; Theoretical: **>1 Mb (ultra-long reads)**                     |
| **Devices**              | **MinION** (USB, portable, field use) <br> **PromethION** (high-throughput, lab-scale) |
| **Unique features**      | Direct RNA sequencing (no cDNA needed); real-time data streaming                       |
| **Pros**                 | Very low-cost device, portable (field sequencing), ultra-long reads, RNA direct sequencing |
| **Cons**                 | Higher error rates (~5–20%), throughput less stable, library prep still required, not good for degraded samples |

<img width="306" height="363" alt="Screenshot 2025-08-29 154854" src="https://github.com/user-attachments/assets/8031d5c6-b0e8-4b6b-a461-6c44e1694d4e" />

<img width="331" height="385" alt="Screenshot 2025-08-29 154946" src="https://github.com/user-attachments/assets/2839c194-5057-4551-9199-14d52d075f80" />

---

## 3. Raw data processing and file formats

**NGS** produces raw signal signal data that must be converted into reads with per-base quality scores.
Step include:
*-1 Basecalling*: conversion of raw signals into nucleotides
*-2 Read filtering*: remove poor-quality sequences
    Not all reads are equally useful. Filtering is so crucial: remove polyclonal reads, discarded reads<25 bases too short for reliable alignment, trim low-quality regions     (sliding windows/moving windows) and adjust parameters for low-complexity libraries.
*-3 Alignment*: mapping reads to a reference genome: once sequencing data (FASTQ) are aligned against a reference genome teh results are stored in SAM and BAM files.
*-4 Variant Calling*: identifying SNPs, indels or structural variants: After sequencing, alignment, QC, and filtering, the next step is variant calling. The goal is to         identify genetic differences (SNPs, indels, SVs) compared to a reference genome. Main software: GATK (Genome Analysis Toolkit) 

<img width="365" height="329" alt="Screenshot 2025-08-29 155044" src="https://github.com/user-attachments/assets/2895db94-2ee0-4a1e-bddf-cede451780c0" />

---

* **FASTA** format: plain text with nucleotide sequence (ACGT..) + identifier. No quality information.
<img width="374" height="334" alt="Screenshot 2025-08-29 154539" src="https://github.com/user-attachments/assets/052e8b94-4f80-4cf7-9771-c0b107412b9d" />

  
* **FASTQ**: raw reads + quality-->FASTA + Phred quality scores per base
    Four lines per entry: Sequence ID (header); Nucleotide sequence; Separator (+); Quality scores (ASCII symbols)
  **quality scores (Phred scale)**: Q=-10log10(e) where Q stay for quality score and e for probability of incorrect base call.
  ex. Q20-->1 error in 100 bases (99% accuracy); Q30-->1 error in 1000 bases (99.9% accuracy)
  Ion torrent error rates Q20 due to homopolymers improvinf towards Q30.
  FASTQ-->unaligned reads.
<img width="383" height="333" alt="Screenshot 2025-08-29 154602" src="https://github.com/user-attachments/assets/58275b5d-feb8-43a9-bb9c-998422234d32" />


  **Quality control (QC)** is a fondamental step in NGS workflows: ensures that the raw data (FASTQ...) are reliable; identifies biases introduced by library preparation,    sequencing platform or sample quality; prevent False positives in downstream variant discovery.
  Why matters: poor-quality reads, GC biases misinterpretation of gene regions, duplicates..
  **Tool for QC**: FastQC java-based, modular, HTML report, and Prinseq efficient for trimming and filtering.

  **FASTQC** -->provides modular analyses: import data from BAM, SAM, FASTQ, quick overview of potential problems, summary graphs, export results and tables and work         offline.
  -Basic statistic (sequence count, length distribution)
  -Per-Base sequence quality (signal): boxplots per position along the read, aggregated phred scores show accuracy, Q>=20 and green zone good orange warning; quality drops    towards the end of reads
  
  <img width="321" height="359" alt="Screenshot 2025-08-29 155146" src="https://github.com/user-attachments/assets/d51e1291-24b1-492f-bbf9-2c21c4497565" />

  -Per-sequence quality score:in the X Phred in the Y the number of reads, shows distribution of read-level quality, reads with mean Q<20 should be discarded.
  
  <img width="320" height="367" alt="Screenshot 2025-08-29 155154" src="https://github.com/user-attachments/assets/e0b35112-f141-49b5-9d21-997e3ca4619e" />

  -Per-base sequence content: shows frequency of ATCG at each position; stable percentages across read length; fluctuantions at the start/end: sequencing errors or library    prerp bias.

<img width="338" height="323" alt="Screenshot 2025-08-29 155201" src="https://github.com/user-attachments/assets/9889b602-bf80-47ec-88b2-c882c72ea987" />


  -GC content: distribution of GC content per read; compared against theoretical distribution of the target genome; deviations indicate biases in library prep and             contaminants; if multiple peaks-->contaminants or mixed DNA sources.

  <img width="303" height="314" alt="Screenshot 2025-08-29 155206" src="https://github.com/user-attachments/assets/ff8ea1f5-91da-43e0-b190-6953da2551eb" />

  
  -Sequence duplication levels: most reads are unique (expected); high duplication bias in library prep; amplicon sequencing duplication expected; duplicated reads can        bias varaint calls, must be removed; after deduplication most reads should be unique.

  <img width="397" height="340" alt="Screenshot 2025-08-29 155220" src="https://github.com/user-attachments/assets/52ea32a0-555b-4a9f-afb3-749cdb850387" />


  **Trimming**-->removing low-quality bases increases dataset reliability. 2 main approaches:
  -Threshold-based trimming: define Q threshold (Q20), removes bases below threshold until high quality base reached
  -Window-based trimming (preferred): define window size (5 nt), calculate average quality in the window, if average<threshold you identify trim region, retains more nt by    smoothing local fluctuations.

<img width="361" height="301" alt="Screenshot 2025-08-29 155254" src="https://github.com/user-attachments/assets/9a0b4bbc-cfdd-435d-af43-7fe4ace55bc9" />


  *PIPELINE OF QC*-->run FASTQC; apply trimming/filtering; remove low-quality bases and duplicated reads; re-run FastQC; only then proceed with alignment and variant         discovery
  
---

* **SAM/BAM**: alignments, CIGAR strings.
  **BAM** stay for binary alignment map. compressed binary format with both reads + alignment information; requires index file; not human readable. It is the binary          version of SAM.
  **SAM** stay for sequence alignment map: a text based format, human-readable.
  A SAM file has 2 sections:
  **a) Header**
  Lines start with @.
  Contains metadata about alignment:
  Reference sequence IDs & lengths.
  Alignment program & command line (@PG).

  **b) Alignment Section**
  Each read is represented by a row with 11 mandatory fields:
  QNAME → read identifier.
  
  FLAG → integer encoding alignment properties.
         The FLAG encodes properties of each read as integers.
         Examples: 4 → read is unmapped. 256 → secondary alignment. 1024 → PCR duplicate.
         Multiple values can be combined.
  
  RNAME → reference sequence name (chromosome).
  POS → starting position on reference.
  MAPQ → mapping quality score.
  
  CIGAR → compact description of alignment.
          CIGAR = Compact Idiosyncratic Gapped Alignment Report. Describes how each read aligns to the reference.
          Common symbols: M → match/mismatch (aligned). I → insertion (present in read, absent in reference). D → deletion (absent in read, present in reference).
                          S → soft clipping (part of the read not aligned).
  
  RNEXT → mate/next read info.
  PNEXT → position of mate/next read.
  TLEN → observed template length.
  SEQ → read sequence.
  QUAL → base quality (Phred).
  
<img width="398" height="347" alt="Screenshot 2025-08-29 155309" src="https://github.com/user-attachments/assets/38c9c59f-2460-495e-85d7-7a8333c453b6" />

  
  BAM and SAM are standard formats in genomics and widely used in downstream analysis.

  *PIPELINE*:start form FASTQ file-->align read to reference genome with BWA-MEM-->the output is the BAM file (compressed alignment)-->BAM can be converted back to SAM.
  
  **Quality control** continues after alignment:
  -Mapping Quality (MAPQ)-->Probability that read is correctly mapped.
   MAPQ = 0 → read could map to multiple locations (e.g., repetitive elements, assembly errors).
  Reads with low MAPQ often discarded.
  -Duplicate Removal-->Duplicates come from PCR amplification.
  They don’t add information, but inflate coverage.
  Can create false positives in variant calling.
  Tools: Picard for duplicate marking/removal.
  PCR-free library preparation avoids this issue (requires more input DNA).
  
---

**Calling a variant** requires careful evaluation of multiple factors: Base call quality of each supporting base. Mapping quality (MQ) of aligned reads. Sequencing depth → minimum number of reads supporting the variant. Proximity to indels or homopolymer runs (error-prone regions).
Single-sample vs multi-sample calling. More samples = higher confidence.

**variant calling strategies**: 
a) Per-sample calling
   Each sample analyzed separately.
   Produces one VCF file per sample.
b) Joint calling--> All samples analyzed together. Produces one merged VCF file.
   Advantages:
   Increases statistical power.
   Better detection of heterozygous variants.
   Reduces false negatives from low coverage.
   Recommended: joint calling for population studies.
  
* **VCF**: SNPs, indels, SVs.
  VCF stay for variant call format: text format for variants (SNPs, indels, structural differences)
  A VCF file has two sections:
  a) Header
  Lines start with ##.
  Define metadata, filters, INFO fields, FORMAT descriptors.
  Last header line starts with #CHROM → column names.
  
  b) Variant Records
  Each row represents one variant.
  Columns include:
  CHROM → Chromosome ID.
  POS → Position of variant.
  ID → Variant ID (e.g., dbSNP rs number, or . if novel).
  REF → Reference allele.
  ALT → Alternative allele(s).
  QUAL → Phred-scaled variant quality score.
  FILTER → PASS/FAIL (e.g., “PASS”, “q10”).
  INFO → Additional annotations (DP = depth, AF = allele frequency, etc.).
  FORMAT → Defines per-sample genotype fields (e.g., GT, DP, GQ).
  10+. Samples → Genotypes and metrics for each sample.
  Genotype Codes
  0|0 → homozygous reference.
  0|1 or 1|0 → heterozygous.
  1|1 → homozygous alternate.
  If multiple alternate alleles exist: 2, 3, etc.
  
<img width="306" height="374" alt="Screenshot 2025-08-29 155325" src="https://github.com/user-attachments/assets/bf17d08a-9184-4c0f-9592-e1ae3e69e8ce" />


  Visualization: **IGV (Integrative Genomics Viewer)** → gold standard for inspecting alignments.
  Inputs: reference genome (FASTA), alignment (BAM), and variants (VCF).
  Features:
  Visualize reads aligned at a given locus.
  Detect true variants vs sequencing errors.
  Depth bar shows number of supporting reads.
  Color-coding distinguishes reference vs alternate alleles.

  Other Tools for Variant Annotation:
  **Ensembl Genome Browser**-->Allows visualization of genomes, genes, transcripts, proteins, and variants.
  Example: searching the KMO gene gives: Chromosomal location. Gene description. Transcript IDs and number. Encoded protein + UniProt link. Variant table for known           polymorphisms (dbSNP, Ensembl).

  **Variant Effect Predictor (VEP)**-->Provided by Ensembl.
  Takes VCF as input.
  Outputs: Genomic location of variants. Variant type (missense, synonymous, frameshift…). Functional consequences on gene/protein.
  **dbSNP**-->Repository of known SNPs and small indels. Variants deposited with IDs (rsXXXX). Can be cross-referenced to validate if a variant is novel.


  **Types of variants**-->
  
  <img width="385" height="352" alt="Screenshot 2025-08-29 155939" src="https://github.com/user-attachments/assets/7b0669a8-c647-48ee-99b6-cd02bf71593b" />
  
<img width="365" height="371" alt="Screenshot 2025-08-29 155021" src="https://github.com/user-attachments/assets/4463cba2-bda9-451d-937f-8155b868e7e1" />


  **Small Variants**: SNPs (Single Nucleotide Polymorphisms), Substitution of one base with another.Must occur in ≥1% of the population → polymorphism.
  If <1%, it is considered a mutation.
  Indels: Small insertions or deletions.
  Substitutions: Multiple bases replaced by different ones.

  **Structural Variants**--> CNVs (Copy Number Variants) → duplicated or deleted regions.
  Inversions → DNA region flipped in orientation.
  Translocations → block of DNA moved to another chromosome or region.

---
  
* **GFF/GTF, BED**: annotations & intervals.
  BED stay for browser extensible data: tab-delimited text format defining genomic features.
  **GFF3/GTF** = rich annotation format (genes, transcripts, exons, CDS, regulatory elements).
  -->Purpose: Stores annotations (genes, exons, CDS, regulatory features, etc.) on a reference genome.
     Format: 9 mandatory tab-delimited columns
  seqid → chromosome or scaffold name
  source → annotation source (e.g. Ensembl, maker, augustus)
  type → feature type (gene, exon, CDS, mRNA, repeat…)
  start → start coordinate (1-based)
  end → end coordinate
  score → numerical value (or “.” if not used)
  strand → + or – (strand of feature)
  phase → for CDS: 0, 1, or 2 (frame of translation start)
  attributes → semicolon-separated key=value pairs (e.g., ID=gene1;Name=BRCA1)

<img width="378" height="381" alt="Screenshot 2025-08-29 155659" src="https://github.com/user-attachments/assets/1c7ecfe0-995b-41f6-9199-da99f96b7ac9" />

  **BED** = simple intervals (coordinates, lightweight for browsers like UCSC/IGV).
  -->Purpose: Defines genomic intervals (regions of interest).
  Format: Very light, 3 mandatory columns, max 12 optional.
  Optional (up to 12 fields): name, score, strand, thickStart, thickEnd, itemRgb, blockCount, etc.
  
<img width="348" height="323" alt="Screenshot 2025-08-29 155707" src="https://github.com/user-attachments/assets/e12cb301-d0be-4674-809d-c70aad884865" />


Raw sequencing output can be huge (terabytes); filtering and compression reduce size dramatically (2.5 TB raw data-->30 GB FASTQ); effcient storage and file format choice are critical for downstream analysis.

| Format      | Purpose / Content                                     | Structure / Key Fields                                                                 | Example (simplified)                                   |
|-------------|-------------------------------------------------------|----------------------------------------------------------------------------------------|--------------------------------------------------------|
| **FASTA**   | Stores raw nucleotide or protein sequences (no quality) | `>identifier` + sequence lines                                                          | `>seq1` <br> `ATGCCGTA...`                             |
| **FASTQ**   | Raw reads + **per-base quality scores**               | 4 lines/entry: <br> 1. `@ID` <br> 2. sequence <br> 3. `+` <br> 4. quality (ASCII, Phred) | `@read1` <br> `ATGCC` <br> `+` <br> `IIIII`            |
| **SAM**     | Text format: aligned reads to reference genome        | Header (`@SQ`, `@PG` …) + alignment section (11 mandatory fields: QNAME, FLAG, RNAME, POS, MAPQ, CIGAR, RNEXT, PNEXT, TLEN, SEQ, QUAL) | `read1  0  chr1  100  60  50M  *  0  0  ATGCC...  IIIII` |
| **BAM**     | Binary compressed version of SAM                     | Same fields as SAM but **binary (indexed)** for speed & storage                          | (binary, not human-readable)                          |
| **VCF**     | Variant calls: SNPs, indels, SVs                     | Header (`##` metadata, `#CHROM` cols) + records (CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT, samples) | `chr1  123456  rs123  A  G  99  PASS  DP=30  GT:DP  0/1:15` |
| **GFF3**    | Genome annotations: genes, exons, CDS, regulatory features | 9 tab-delimited cols: seqid, source, type, start, end, score, strand, phase, attributes | `chr1 Ensembl gene 1000 5000 . + . ID=gene1;Name=BRCA1` |
| **BED**     | Genomic intervals (lightweight)                      | 3 required cols: chrom, start, end (up to 12 optional: name, score, strand, color, blocks) | `chr1  1000  5000  BRCA1`                              |

---

**Before sequencing: DNA Quality Assessment & Sequencing Coverage**-->

1. Quality Assessment of DNA
Before starting NGS library preparation, DNA quality must be carefully evaluated.
DNA purity and integrity affect downstream steps such as library prep, amplification, and sequencing efficiency.

**Purity Ratios (Spectrophotometry)**
A260/280 ratio
Expected: ~1.8
High values → RNA contamination
Low values → protein contamination, residual phenol, or very low DNA concentration (<10 ng/µl)
Influenced by pH (acidic solution ↓ ratio; basic solution ↑ ratio)
A260/230 ratio
Expected: 2.0 – 2.2
Low values → carbohydrate carryover (plants), residual phenol, guanidine, glycogen
High values → improper blank solution

**DNA Integrity**
Checked by agarose gel electrophoresis:
High-quality DNA → sharp band, high molecular weight (upper gel region)
Degraded DNA → smeared band or lower fragments
Low integrity does not prevent sequencing, but long-read sequencing (Nanopore, PacBio) requires intact long fragments.

2. Sample Requirements for NGS
Companies providing sequencing services usually require:
Condition: genomic DNA
Quantity:
≥ 1 µg (small fragment library)
1.5 µg (PCR-free library)
Concentration: ≥ 12.5 ng/µl
Purity: A260/280 between 1.8 – 2.0
If requirements fail, sequencing may still be possible, but with increased risk of poor results.

3. Sequencing Coverage
**Coverage (breadth)** → % of target bases sequenced at least “X” times.
Depth of coverage → average number of times each base is sequenced.
Depth of Coverage} = LN/G

Where:
L = read length
N = number of reads
G = haploid genome length
Expressed as “X-fold” (e.g., 10x, 20x, 40x).
Examples
5x coverage → each base sequenced ~5 times (low reliability, missing variants).
20x coverage → sufficient for variant discovery (SNPs/indels).
40x coverage → high-quality genome resequencing (each base seen ≥40 times).
Genome assembly → requires much higher depth and often long-read technologies.

4. Why Coverage Matters
Variant discovery → ≥20x needed for reliable SNP/indel detection.
Clinical genomics → high coverage required to avoid missing pathogenic variants.
De novo assembly → high depth + long reads needed to resolve repeats and structural variants.

---

## 4. Genome Assembly

Genome assembly is the process of reconstructing the complete genomic sequence of an organism from millions of sequencing reads. Modern sequencing platforms produce vast numbers of short or long fragments, which must be computationally pieced together to rebuild the genome.

After sequencing, genome assembly is required when no reference genome is available to perform alignments.

Assembly can serve multiple purposes:

Create a reference genome for a species that does not yet have one.
Improve an existing reference genome.
Genome assembly is strictly related to genome annotation, since it is fundamental to know where reads belong in order to reconstruct the genome.
Because de novo genome assembly requires high computational power, time and expert teams, researchers must evaluate in advance whether it is truly needed or if an existing reference genome can be used.

### Shotgun Sequencing

The most common strategy is **shotgun sequencing**, where DNA is randomly fragmented into smaller pieces that are sequenced in parallel. This provides unbiased coverage of the genome but creates the challenge of correctly reassembling overlapping fragments.

DNA is extracted and fragmented to build a library.
Overlaps between fragments are analyzed to reconstruct long stretches of DNA.

*Hierarchical Shotgun Sequencing*

<img width="386" height="355" alt="Screenshot 2025-08-29 155412" src="https://github.com/user-attachments/assets/e03687f2-e993-426c-aef7-ed3161193b37" />

Used in early genome projects.
Genome cut into large fragments (contigs) → inserted into BAC libraries (Bacterial Artificial Chromosomes).
Each BAC is shotgun-sequenced and assembled.
Individual BAC assemblies are aligned to recover the full genome.
Nowadays: replaced by direct whole-genome shotgun sequencing thanks to computational power.

### WORKFLOW
General steps
Sample collection → obtain representative tissues. Some can be saved for RNA-seq (useful for annotation).
DNA extraction → decision between long-read sequencing (better continuity, e.g., PacBio, Nanopore) and short-read sequencing (higher coverage but fragmented, e.g., Illumina).
Genome assembly pipeline → composed of 10 steps (see below)-->

-Gather information about the target genome
Genome size, heterozygosity, repeats, ploidy.
Methods to estimate size:
Flow cytometry → provides C-value.
Databases → use closely related species.
K-mer frequency distribution → counts subsequences in reads.

-Extract high-quality DNA
Long stretches (up to 20 kbp) for long-read sequencing.
Short reads used for degraded/ancient DNA.
Remove contaminants (mitochondrial, chloroplast, microbial DNA).

-Design experimental workflow
Define experimental goals and sequencing depth (≥ 60x recommended).
Choose library type (single-end, paired-end, mate-pair, PCR-free).
Consider costs (PacBio vs Illumina).
Choose between de novo and reference-guided assembly.

-Choose sequencing technology & library preparation
Illumina: cheap, fast, high-throughput short reads.
PacBio/Nanopore: expensive but long reads.
Hybrid approach: correct long-read errors with short reads.
Use additional technologies (Hi-C, BioNano) for scaffolding.

-Ensure computational resources
Example: diploid genome (~1 Gb) → 96 CPUs, 1 TB RAM, 3 TB local storage, 10 TB shared.
Cloud services (e.g., AWS) often required.

-Assemble the genome
Preprocess reads (QC, filtering).
Use algorithms (Greedy, OLC, De Bruijn graph, String graph, Hybrid).

-Assembly polishing
Correct sequencing errors, improve consensus sequence.
Errors can come from repeats, misoriented contigs, or sequencing artefacts.

-Check assembly quality
Metrics: assembly size, N50, completeness.
N50 = contig length at which 50% of the assembly is contained in contigs of that size or longer.

-Scaffolding & gap filling
Merge contigs using paired-end or mate-pair reads.
Use genetic linkage maps or reference-guided scaffolding.

-Genome annotation
Map genes and regulatory elements using RNA-seq and other functional data.

### Assembly Algorithms
Different computational strategies are used depending on read length and sequencing technology:

| Feature                | **Overlap–Layout–Consensus (OLC)**                                | **de Bruijn Graph**                                       |
|-------------------------|-------------------------------------------------------------------|-----------------------------------------------------------|
| Input                   | Long reads (Sanger, PacBio, Nanopore)                            | Short reads (Illumina)                                    |
| Principle               | Reads compared pairwise for overlaps → graph of overlaps → consensus | Reads decomposed into **k-mers**; overlaps represented as graph edges |
| Pros                    | Accurate with long reads; good for smaller datasets              | Efficient with huge short-read datasets; scales well       |
| Cons                    | Computationally expensive (all-vs-all overlaps); not efficient for short reads | Struggles with sequencing errors, repeat resolution depends on *k* size |

*Greedy algorithms* → iterative merging of overlapping reads; too simple for NGS.
Graph theory is introduced to solve complexity:
Vertices = reads/k-mers, Edges = overlaps.
Types of graphs: undirected, directed, multigraph.
In-degree / out-degree concepts used in Eulerian paths.

*OLC*--> Compute overlaps between reads.
Build graph and create contigs.
Infer consensus sequence.
Accurate but computationally expensive for large NGS datasets.

<img width="367" height="301" alt="Screenshot 2025-08-29 155505" src="https://github.com/user-attachments/assets/3312fa30-ae2f-44cf-9abb-3f77d7cb8c9c" />


*De Bruijn Graph (DBG)*-->
Reads decomposed into k-mers.
Vertices = (k–1)-mers; Edges = k-mers.
Efficient graph traversal reconstructs contigs.
Errors can cause unique erroneous k-mers; choice of k is crucial.
Uses Eulerian paths (each edge visited once).
Repeats lead to branched structures, resolved by longer libraries or paired/mate-pair reads.

<img width="332" height="380" alt="Screenshot 2025-08-29 155518" src="https://github.com/user-attachments/assets/38d0167f-7184-419d-bda7-1fb3f5d0fbe9" />


### k-mers and Genome Size Estimation
**k-mers** (substrings of length *k*) are the foundation of short-read assembly. Counting the frequency of all k-mers across sequencing reads produces a distribution that can be used to estimate:
- **Genome size**:  
 Genome size = Total k-mers/Peak depth
- **Heterozygosity**: heterozygous genomes show two peaks (diploid k-mer distribution).  
- **Repeat content**: repetitive sequences distort the k-mer curve and increase multiplicity.
  
  <img width="396" height="358" alt="Screenshot 2025-08-29 155438" src="https://github.com/user-attachments/assets/3aa32a8d-933d-49bb-be43-b1af6b9bdb27" />


### C-Value (Absolute Genome Size)
The **C-value** is the amount of DNA contained in a haploid nucleus (1C, e.g., in gametes). It represents the absolute genome size, measured in base pairs or picograms of DNA.  
- Example: the human haploid genome has a **C-value ≈ 0.978 × 10⁹ bp (~3.2 pg DNA)**.  
- C-value is species-specific and important in planning sequencing projects (coverage and depth requirements).
  
<img width="377" height="325" alt="Screenshot 2025-08-29 155428" src="https://github.com/user-attachments/assets/d48c2053-41f4-438c-bf0d-e950a434fc75" />

### Scaffolding
Once contigs (continuous assembled sequences) are generated, **scaffolding** uses additional information to connect and order them:
- **Paired-end reads**: link short fragments separated by a few hundred bp  
- **Mate-pair libraries**: provide longer inserts (2–20 kb) for spanning repeats  
- **Long reads (PacBio, Nanopore)**: bridge large gaps and improve contiguity  


### Assembly Quality Metrics

| Metric       | Definition & Details                                                                                                      |
|--------------|---------------------------------------------------------------------------------------------------------------------------|
| **N50**      | A contiguity statistic. All contigs/scaffolds are ordered from longest to shortest, and lengths are cumulatively summed until 50% of the total assembly size is reached. The length of the contig at this point is the **N50**. <br> Example: if genome = 3 Gb, and adding contigs from longest to shortest passes 1.5 Gb at contig length = 5 Mb, then **N50 = 5 Mb**. |
| **Coverage** | Average sequencing depth per base (e.g. 20×, 40×). Calculated as: <br> \[ \text{Coverage} = \frac{L \times N}{G} \] where *L* = read length, *N* = number of reads, *G* = genome length. Higher coverage improves accuracy and completeness. |
| **BUSCO**    | **Benchmarking Universal Single-Copy Orthologs**: searches for a set of evolutionarily conserved orthologous genes expected in a lineage. Results are reported as percentages of: <br> - **Complete (C)** → full-length ortholog found <br> - **Single-copy (S)** → present once <br> - **Duplicated (D)** → multiple copies detected <br> - **Fragmented (F)** → partial genes <br> - **Missing (M)** → not found <br> Example: BUSCO = 95% Complete (90% S, 5% D), 3% Fragmented, 2% Missing → indicates a very complete assembly. |

<img width="369" height="315" alt="Screenshot 2025-08-29 155545" src="https://github.com/user-attachments/assets/e3ce6a26-b72c-4f7a-adcb-5de1b29101cb" />


**BUSCO (Benchmarking Universal Single-Copy Orthologs)**:
Evaluates genome completeness by searching for conserved single-copy genes.
If genes appear more than once → possible misassembly.
Useful for both assembly and annotation quality control.
**In summary**, genome assembly integrates sequencing technologies, graph-based algorithms, k-mer analysis for genome size and C-value estimation, scaffolding strategies, and rigorous quality metrics (N50, coverage, BUSCO) to transform fragmented reads into a biologically meaningful genome sequence.

<img width="326" height="382" alt="Screenshot 2025-08-29 155616" src="https://github.com/user-attachments/assets/656db6e1-3d88-4a33-bce4-244b9263d634" />

---
## GALAXY
-->workflow management system, accesible, reproducible, transparent

<img width="312" height="381" alt="Screenshot 2025-08-29 155347" src="https://github.com/user-attachments/assets/450c40a8-ea4f-4c4a-bf5c-d7f1f34278af" />

---

## 5. Genome Annotation

<img width="344" height="325" alt="Screenshot 2025-08-29 155640" src="https://github.com/user-attachments/assets/60e2878d-d5c1-4089-b854-8c31f99496d3" />

Genome annotation is the process of adding biological meaning to an assembled genome.
It involves:

Identifying repeats.
Locating genes (structural annotation).
Assigning functions (functional annotation).
Using intrinsic, extrinsic, or combined prediction methods.
Requiring manual curation to achieve gold-standard annotations.
Submitting final data in standardized formats to public repositories.

Once a genome has been assembled, the next step is **annotation**, the process of identifying and describing the functional elements within the sequence. Annotation transforms a raw collection of contigs and scaffolds into a biologically meaningful map of genes, repeats, regulatory regions, and other features.  

*Genome Annotation Steps-->*
Annotation = assigning biological roles to genomic sequences.
Two main phases:

Feature identification (using algorithms and pipelines).
Annotation (synthesizing data).
Main annotation workflow
Annotation of repetitive elements (REs).
Structural annotation → locating genes.
Functional annotation → assigning functions to genes.
Data submission, maintenance, updates.

---

### Prokaryote vs Eukaryote Annotation
Prokaryotes
Genes = long ORFs (>50 codons).
Low intergenic DNA.
No introns.
Easy to detect!
Eukaryotes
Genes interrupted by introns.
High intergenic DNA (e.g., 62% in humans).
Complex regulatory features and UTRs.
Need to consider codon bias, intron–exon boundaries, regulatory motifs.

### Transcript Information
Short-read RNA-seq: high throughput but fragmented → harder to map accurately.
Long-read RNA-seq (PacBio/Nanopore): captures full-length transcripts, resolves isoforms.

### Annotation Targets
Protein-coding genes (with ORFs).
Non-coding RNAs (lncRNAs, snoRNAs, miRNAs).
Pseudogenes → annotated using paralog/ortholog homology (tools: PseudoPipe).

---
### Repeat Annotation
A large proportion of most eukaryotic genomes is composed of **repetitive elements**, including transposable elements (TEs), tandem repeats, and low-complexity regions. Detecting and masking repeats is crucial because they can cause false gene predictions and complicate downstream analyses.  

- **RepeatMasker** → the most widely used tool, screens DNA sequences for interspersed repeats and low complexity sequences using curated repeat libraries (e.g., RepBase, Dfam).  
- **Dfam** → curated database of transposable element families, based on profile HMMs.  
- **TEannot (REPET pipeline)** → specialized in de novo transposable element discovery and annotation.  

Repeat annotation ensures that repetitive sequences are catalogued, masked when necessary, and separated from protein-coding genes. This prevents false gene calls and improves downstream prediction.  

### REs
Types
Low-complexity sequences: e.g., homopolymeric runs (AAAAA).
Transposable elements (TEs):
Class I: RNA intermediate, copy & paste.
Class II: DNA intermediate, cut & paste.

-Challenges
Borders are not always clear.
Nested repeats (repeats inside repeats).
Huge variety of classes → time-consuming.

-Tools and Approaches
Library-based (homology)
Compare genome vs curated libraries (species-specific or general).
Examples: REPBASE, RepeatMasker (uses BLAST, HMMs like NHMMer).
De novo discovery
Identify candidate repeats based on structural similarity.
Mixed approaches

-Databases
Dfam:
Open database of TE profiles (HMMs + consensus sequences).
Built from multiple sequence alignments → seed alignments → profile HMMs.
REPET package:
Tedenovo → detects TEs.
Teannot → annotates TEs.
Masking repeats
Replace nucleotides with N (hard masking) or lowercase letters (soft masking).
Importance: prevents spurious BLAST hits and false gene annotations.

**Table – Repeat Annotation Tools**

| Tool/DB        | Principle                              | Output                           | Applications                       |
|----------------|----------------------------------------|----------------------------------|------------------------------------|
| **RepeatMasker** | Uses libraries (RepBase, Dfam) to find interspersed repeats | Annotated repeat-masked genome   | Genome masking, TE annotation      |
| **Dfam**       | Profile HMMs for TE families           | TE classification and consensus  | Reference for repeat detection     |
| **TEannot**    | De novo discovery + homology           | Custom TE library + annotation   | Non-model species, new TE families |

---

### Gene Prediction and Gene Models
The central task of genome annotation is to define **gene models**, which describe the structure of genes, including exons, introns, untranslated regions (UTRs), and coding sequences (CDS).  

Several strategies exist, often combined for best accuracy:  

1. **Ab initio prediction** → uses intrinsic sequence signals (start/stop codons, splice sites, codon bias).  
   - **AUGUSTUS** is the gold-standard tool; trained on reference genomes, it predicts complete gene structures from sequence alone.  
2. **Homology-based prediction** → aligns known proteins/transcripts from related species to infer gene models. Tools: **BLAST, Exonerate**.  
3. **Transcript evidence** → integrates RNA-seq data; mapped reads reveal expressed exons, introns, and splice junctions.  
4. **Integrative approaches** → pipelines like **MAKER, BRAKER2** combine ab initio predictions with transcript and homology evidence, creating consensus models.  

This **multi-layered approach** produces more reliable gene models, capturing both coding and non-coding RNAs.  

**Table – Ab initio vs Extrinsic/Integrative Approaches**

| Approach / Tool   | Input Data                        | Strengths                                            | Limitations                          |
|-------------------|-----------------------------------|------------------------------------------------------|--------------------------------------|
| **Ab initio (AUGUSTUS)** | DNA sequence only                  | Detects novel genes, independent of external data     | False positives in repetitive regions, needs training |
| **Homology-based** | Known proteins/transcripts        | High accuracy for conserved genes, useful in non-models | Misses species-specific/novel genes   |
| **Transcript evidence** | RNA-seq reads, ESTs               | Captures real expression, splicing isoforms           | Limited to expressed genes, condition-dependent |
| **Integrative (MAKER, BRAKER2)** | Combines ab initio + homology + RNA-seq | Highest reliability, consensus models, widely used    | Computationally intensive, needs multiple data types |

*Intrinsic (Ab initio)*
Relies only on genome sequence itself.
Uses mathematical models trained per genome.
Predicts ORFs (start codon → stop codon).
Pros: high sensitivity with enough training data.
Cons: lower accuracy for intron–exon boundaries.

*Extrinsic*
Relies on databases of transcripts or proteins.
Uses homology to known genes.
Sources: RefSeq, UniProt, NCBI nr, RNA-seq data.
Pros: universally applicable.
Cons: risk of missing novel, species-specific genes.

*Combiners*
Integrate ab initio + extrinsic approaches.
Examples: AUGUSTUS (can work in both ways).
Allow one type of evidence to override another if it improves prediction.

---

### Functional Annotation
After structural annotation, the next step is to assign **biological meaning** to genes and transcripts. This is achieved through cross-referencing with biological databases:  

- **Gene Ontology (GO)** → assigns terms for *Molecular Function, Biological Process, Cellular Component*.  
- **KEGG (Kyoto Encyclopedia of Genes and Genomes)** → maps genes into pathways (metabolism, signaling).  
- **Protein domains (Pfam, InterPro)** → detect conserved motifs, functional domains, and evolutionary relationships.  

Functional annotation enables downstream interpretation, linking raw sequence to biology, pathways, and phenotypes.  

---

### Genome Annotation – Pipeline Overview

1. **Repeat Annotation** → identify and mask repeats (RepeatMasker, Dfam, TEannot).  
2. **Structural Annotation** → predict gene models (ab initio, homology, transcript evidence, integrative pipelines).  
3. **Functional Annotation** → assign GO terms, KEGG pathways, Pfam/InterPro domains.  
4. **Curation** → manual review of key genes, correction of misannotations.  

---

### In Summary
Genome annotation transforms raw assemblies into functional blueprints of organisms. It requires **repeat masking**, **robust gene prediction**, and **functional annotation**. High-quality annotation is critical for:  
- **Comparative genomics** (ortholog/paralog detection, synteny).  
- **Transcriptomics** (RNA-seq alignment and quantification).  
- **Functional studies** (gene discovery, pathway analysis).  
- **Applied biotechnology** (breeding, engineering, synthetic biology).  

---

## 6. Specialized Sequencing

In addition to whole-genome sequencing, a range of **specialized approaches** targets specific biological questions. These methods focus on defined regions or molecular layers and complement standard WGS/RNA-seq.

---

### aCGH (Array Comparative Genomic Hybridization)

<img width="359" height="328" alt="Screenshot 2025-08-29 160018" src="https://github.com/user-attachments/assets/6a4cb9a6-0961-4584-a678-cccd3a93b375" />

**Principle**: Detects **Copy Number Variations (CNVs)** using hybridization on a microarray. Test and reference DNA are labeled with different fluorescent dyes, co-hybridized to probes on a chip, and the fluorescence ratio indicates gains or losses.  

**Pipeline**:
1. **DNA extraction** (test + reference).  
2. **Labeling** → test (red) and control (green).  
3. **Co-hybridization** on **microarray chip** with thousands of probes.  
4. **Washing & scanning** → measure fluorescence intensity.  
5. **Signal analysis** → log₂ ratio (test/control).  

**CNV signal interpretation**:
- **Ratio ~0** → normal copy number.  
- **Positive log₂ ratio** → duplication/gain.  
- **Negative log₂ ratio** → deletion/loss.  

**Microarray “chip pipeline”**:  
DNA/RNA → labeling → hybridization on chip probes → scanner reads fluorescent signals → computational normalization → intensity plots → interpretation.  

**Table – aCGH Overview**

| Feature        | Details                                                                 |
|----------------|-------------------------------------------------------------------------|
| **Principle**  | Test vs reference DNA → labeled → hybridized on microarray → intensity ratio = CNV |
| **Resolution** | Typically 20–100 kb (depends on probe density)                          |
| **Applications** | CNV detection in cancer genomics, constitutional disorders, prenatal screening |
| **Strengths**  | Genome-wide CNV detection, relatively cheap, robust technology          |
| **Limitations**| Cannot detect balanced rearrangements (inversions/translocations); limited resolution compared to NGS |


### Microarrays and CGH
- **Microarrays**: DNA probes on a chip; hybridization-based detection of gene expression, SNPs, and known variants.
- **Classical CGH (cCGH)**: co-hybridization of tumor vs reference DNA on metaphase chromosomes; low resolution (~5–10 Mb).
- **Array CGH (aCGH)**: CGH performed on microarrays; higher resolution (kb–Mb) for detecting copy number variations (CNVs).

---

### ChIP-seq Pipeline (Chromatin Immunoprecipitation Sequencing)

<img width="172" height="327" alt="Screenshot 2025-08-29 155824" src="https://github.com/user-attachments/assets/f23b81e5-4169-4ab5-bb79-e7188c181957" />

1. **Crosslinking**  
   - Cells treated with formaldehyde → DNA and bound proteins (e.g., transcription factors, histones) are covalently crosslinked.  

2. **Chromatin fragmentation**  
   - DNA sheared into small fragments (100–500 bp) by sonication or enzymatic digestion.  

3. **Immunoprecipitation (IP)**  
   - Antibody specific to protein of interest (TF, histone mark) pulls down DNA-protein complexes.  

4. **Reverse crosslinking & purification**  
   - Crosslinks removed → recover pure DNA fragments bound by protein.  

5. **Library preparation & sequencing**  
   - Recovered DNA fragments ligated with adapters, amplified, and sequenced (usually Illumina).  

6. **Bioinformatics analysis**  
   - Map reads to reference genome (BWA, Bowtie).  
   - Identify **peaks** = enriched regions of protein binding.  
   - Downstream: motif discovery, regulatory network analysis.

---

### Pool-seq

<img width="365" height="390" alt="Screenshot 2025-08-29 155741" src="https://github.com/user-attachments/assets/570b0c55-0e70-483b-9f88-3c9848aa61f1" />


**Principle**: DNA from many individuals is pooled and sequenced together to estimate **allele frequencies**.  

- **Equimolar DNA pool**: each individual contributes the same DNA amount to avoid bias.  
- Provides allele frequency spectra, selective sweeps, and F_ST between populations, but no individual genotypes.  

**Table – Pool-seq Overview**

| Feature        | Details                                                                 |
|----------------|-------------------------------------------------------------------------|
| **Principle**  | Pool DNA from multiple individuals → sequence together                  |
| **Output**     | Allele frequency spectra, population-level diversity                    |
| **Applications** | Population genomics, selective sweep detection, allele frequency estimation |
| **Strengths**  | Cost-effective, fast, scalable                                          |
| **Limitations**| No individual genotypes, sensitive to unequal DNA contributions         |


## Pool-Seq: Sequencing Pools of Individuals

### Concept
- **Pool-seq** = Whole genome sequencing of pools of individuals.  
- Provides a **cost-effective** alternative to sequencing each individual separately.  
- Allows estimation of **allele frequencies** in a population at lower cost.  

### Workflow
1. Extract DNA from each sample.  
2. Mix DNA in **equimolar quantities** (same number of genome copies per individual).  
   - Quantification done with techniques such as spectrometry.  
3. Sequence the pooled DNA.  
4. Map reads to reference genome.  
5. Call SNPs and calculate allele frequencies.  

Limitation: origin of each read (individual) is lost.  
Advantage: much cheaper while still providing accurate allele frequency estimates.

Applications
- **Detecting genetic basis of phenotypes**.  
- Example: *Genetic basis for red coloration in birds*.  
  - Compare two pools: red vs yellow canaries.  
  - Estimate allele frequencies at each variant.  
  - Identify genomic regions with large allele frequency differences → candidate genes.  
- Useful for **extreme phenotype comparison** (e.g., healthy vs diseased populations).  

### Cost Example
- Study: 56 Gbp sequenced for ~$400 (using a non-Illumina machine; Illumina would cost ~$700).  
- Shows strong cost-effectiveness for population-level studies.

## Targeted Sequencing

### Concept
- **Targeted sequencing** = focus on specific genes or regions of interest.  
- More rapid and cost-effective than sequencing entire genomes.  
- Uses **PCR amplification** or **hybridization-based capture** methods.  

### Approaches
- **PCR amplification** + Sanger sequencing → small number of regions.  
- **Ion AmpliSeq™ panels** → target hundreds of genes in one day (Ion PGM system).  
- **Ion TargetSeq™ Enrichment System** → up to ~60 Mb targeted regions (customizable).  

### Applications
- Cancer genomics.  
- Pharmacogenomics.  
- Discovery of novel variants in selected loci.  

## Amplicon Sequencing (Illumina)

- **Amplicon sequencing** = highly targeted, PCR-based enrichment of specific regions.  
- Produces **ultra-deep sequencing of PCR products (amplicons)**.  
- Benefits:
  - Detect **rare somatic mutations** (e.g., in tumor samples).  
  - Analyze **specific bacterial markers** (e.g., 16S rRNA gene for microbial taxonomy).  
- Workflow:
  - Design primers/probes for regions of interest.  
  - PCR amplification.  
  - NGS sequencing of amplicons.  

## Summary
- **Pool-Seq**:  
  - Effective for allele frequency estimation in populations.  
  - Best for studies comparing **extreme phenotypes**.  
  - Low-cost alternative to sequencing individuals separately.  

- **Targeted Sequencing**:  
  - Focus on specific regions/genes of interest.  
  - Cost-effective, customizable (PCR panels, hybrid capture).  
  - Used in cancer genomics, pharmacogenomics, metagenomics.  

- **Amplicon Sequencing**:  
  - Ultra-deep sequencing of PCR products.  
  - Powerful for rare variant detection and microbial community profiling.  
---

### Targeted Sequencing

**Principle**: Focuses sequencing on specific regions instead of the whole genome.  

- **Amplicon sequencing**: PCR amplifies specific loci.  
- **Hybrid capture panels**: probes enrich predefined genomic regions (e.g., cancer or exome panels).  

**Table – Targeted Sequencing Overview**

| Feature        | Details                                                                 |
|----------------|-------------------------------------------------------------------------|
| **Principle**  | Focused sequencing of selected loci (amplicons or capture panels)        |
| **Applications** | Cancer gene panels, rare disease diagnosis, pharmacogenomics          |
| **Strengths**  | High depth of coverage, cost-efficient, tailored to specific questions  |
| **Limitations**| Misses variants outside target, capture bias, design required in advance |


---

### Exome (Whole-Exome Sequencing — WES)


**Principle**: Captures and sequences protein-coding regions (~1–2% of genome).  

**Pipeline**:
1. DNA fragmentation.  
2. Library prep with adapters.  
3. **Hybrid capture** with biotinylated probes against exons.  
4. Wash & pull-down with streptavidin beads.  
5. PCR amplification.  
6. Sequencing → variant calling.  

**Applications**: Mendelian diagnostics, high-depth coding variant discovery, trio analysis.  

**Table – WES Overview**

| Feature        | Details                                                                 |
|----------------|-------------------------------------------------------------------------|
| **Principle**  | Hybrid capture of coding regions (~20,000 genes, ~1–2% of genome)       |
| **Applications** | Clinical diagnostics, Mendelian disorders, coding SNP/indel discovery |
| **Strengths**  | Cheaper than WGS, deeper coverage, focuses on known disease-causing regions |
| **Limitations**| Misses regulatory/non-coding variants, capture bias, uneven coverage    |

# Whole-Exome Sequencing (WES)

<img width="352" height="354" alt="Screenshot 2025-08-29 155801" src="https://github.com/user-attachments/assets/692074a4-1c12-462e-acbe-1a6b23fe2675" />

## Introduction
- **Whole-exome sequencing (WES)** = NGS method targeting the protein-coding regions of the genome.  
- Human exome = **<2% of the genome**, but contains ~85% of known disease-related variants.  
- Cost-effective alternative to **whole-genome sequencing (WGS)**.  
- Typical data: **4–5 Gb per exome vs ~90 Gb per whole human genome**.  

### Advantages
- Identifies variants across a wide range of applications.  
- Comprehensive coverage of coding regions.  
- Smaller and more manageable data sets → easier and faster analysis.  
- Broad applications: population genetics, rare genetic disease, cancer studies.  


## Hybridization Capture (Target Enrichment)
WES relies on **hybridization capture** to enrich coding regions.

**Workflow**:
1. **Library preparation**  
   - DNA fragmented (mechanical/enzymatic).  
   - Sequencing adapters added.  
   - PCR amplification may be used.  
2. **Hybridization with probes**  
   - Biotinylated probes (“baits”) hybridize to coding regions.  
   - Thousands of probes (up to 500k) target ~20,000 coding regions.  
3. **Capture with streptavidin**  
   - Bait–DNA complexes isolated.  
   - Non-target DNA washed away.  
4. **Sequencing**  
   - Only captured exonic fragments sequenced.  

**Key Point**: Probes are designed to bind specific exonic regions, allowing enrichment of the exome before sequencing.

## Applications
- **Rare variant discovery** → multiple unrelated affected individuals sequenced; shared variants prioritized.  
- **Mendelian disorders** → powerful strategy with modest sample sizes.  
- **De novo mutations** → family-based sequencing (parents unaffected, child affected).  
- **Extreme phenotypes** → sequencing individuals at opposite ends of a quantitative trait (e.g., tall vs short).  

## Variant Filtering Strategy
After WES, thousands of variants are typically identified. Filtering is necessary.

### Discrete Filtering
- **Assumption**: any allele present in a reference “filter set” (databases) is unlikely to be causative.  
- Removes common polymorphisms → focuses on novel/rare variants.  
- Only ~2% of SNVs identified by exome sequencing are novel.  
- Greatly reduces candidate genes to a small, high-priority list.  
- Exceptionally powerful for **rare Mendelian disorders**.  

### Filtering Steps
1. Identify variants shared across affected individuals.  
2. Compare against public databases (e.g., 1000 Genomes, gnomAD).  
3. Exclude common polymorphisms.  
4. Prioritize rare, novel variants with potential disease impact.  

## Comparison of Approaches

| Approach | Target | Coverage per base | Data size | Cost | Notes |
|----------|--------|-------------------|-----------|------|-------|
| **Targeted sequencing** | Specific genes/SNPs | Very high | Very small | Lowest | Focused approach, useful in clinical settings |
| **Whole-exome sequencing (WES)** | All coding regions (~2% of genome) | High | Moderate (~4–5 Gb) | Medium | Captures most disease-causing variants |
| **Whole-genome sequencing (WGS)** | Entire genome | Moderate (spread over 3 Gb) | Large (~90 Gb) | Highest | Comprehensive, unbiased |


## Summary
- WES focuses on **exons**, where most disease-related variants occur.  
- **Hybridization capture** enriches exonic DNA for sequencing.  
- Enables discovery of **rare variants** and candidate genes in Mendelian and complex diseases.  
- **Filtering strategies** (database-based, family-based, extreme phenotype comparisons) refine results.  
- Researchers can choose between **targeted panels, WES, or WGS** depending on cost, data needs, and research goals.

---

### Methyl-seq / Bisulfite Sequencing

<img width="338" height="397" alt="Screenshot 2025-08-29 155817" src="https://github.com/user-attachments/assets/c55e5687-120c-4b99-a56c-5deaf9943cf7" />

**Principle**: Bisulfite converts **unmethylated C → U** (read as T), while methylated C remains unchanged. Sequencing reveals methylation at **base resolution**.  

**Pipeline**:
1. High-quality DNA extraction.  
2. Bisulfite treatment (C → U if unmethylated).  
3. Library preparation (adapters, PCR).  
4. Sequencing.  
5. Alignment to reference genome.  
6. **Methylation calling**: compare C vs T reads → % methylation per CpG.  

**Table – Methyl-seq Overview**

| Feature        | Details                                                                 |
|----------------|-------------------------------------------------------------------------|
| **Principle**  | Sodium bisulfite converts unmethylated C to U (→ T in sequencing); methylated C unchanged |
| **Resolution** | Single base (CpG methylation profiles)                                  |
| **Applications** | Epigenetics, imprinting, cancer methylome, developmental biology     |
| **Strengths**  | High resolution, genome-wide, quantitative                             |
| **Limitations**| DNA damage from bisulfite, incomplete conversion artifacts, requires high input |

# Epigenomics and Methyl-Seq

## Introduction
- **Epigenomics** studies modifications that regulate gene expression without altering the DNA sequence.  
- **DNA methylation** (mainly at cytosines in CpG dinucleotides) plays a key role in gene regulation.  
- The **genome sequence** is the same in all cell types, but methylation patterns differ across tissues, cell types, and developmental stages.  

### CpG Islands
- Regions rich in CpG sites.  
- Often associated with gene promoters.  
- Adjacent regions:  
  - **CpG shores** → within 2 kb from CpG island.  
  - **CpG shelves** → within 2–4 kb.  

**Methylation effect**:  
- Hypermethylation of promoter regions often silences gene expression.  
- Hypomethylation can promote gene activation.  

## Bisulfite Sequencing
**Principle**: treatment of DNA with sodium bisulfite converts cytosines, but **methylated cytosines are protected**.

1. **Unmethylated cytosines (C)** → converted to **uracil (U)** → sequenced as **thymine (T)**.  
2. **Methylated cytosines (5mC)** → remain as **cytosine (C)**.  

### Workflow
1. Extract DNA.  
2. Treat with sodium bisulfite.  
3. Sequence DNA (NGS).  
4. Map reads back to reference genome.  

**Interpretation**:
- If base = **C** → cytosine was methylated.  
- If base = **T** → cytosine was unmethylated.  

## Challenges
- A **C→T change** in sequencing can represent:  
  1. **Bisulfite conversion** (true methylation signal).  
  2. A **genuine variant (SNP)**.  
- Distinguishing between these requires careful experimental design and sometimes complementary sequencing approaches.  

## Library Preparation Strategies
- Different library prep methods have been developed to maximize sensitivity and distinguish:  
  - **Variants (SNPs)** vs **Methylation status**.  
  - Whole-genome vs exome-level methylation.  
- Approaches include:  
  - Whole-genome bisulfite sequencing (**WGBS**).  
  - Reduced representation bisulfite sequencing (**RRBS**).  
  - Targeted bisulfite sequencing (specific loci).  


## Applications of Methyl-Seq
- Studying **tissue-specific regulation**.  
- Identifying **epigenetic biomarkers** in cancer.  
- Understanding developmental processes.  
- Disentangling **genetic variants vs epigenetic modifications**.  

## Summary
- DNA methylation is a key epigenetic mechanism regulating gene expression.  
- **Bisulfite sequencing** is the gold standard for detecting methylated cytosines.  
- Distinguishing between **C→T SNPs** and **C→T due to bisulfite conversion** is critical.  
- Several strategies (WGBS, RRBS, targeted Methyl-Seq) exist to explore methylation patterns across the genome.
  
---

### RNA-seq
<img width="302" height="379" alt="Screenshot 2025-08-29 155854" src="https://github.com/user-attachments/assets/99f5a05b-a41b-45ca-9c6b-5059224a475d" />
<img width="332" height="364" alt="Screenshot 2025-08-29 155907" src="https://github.com/user-attachments/assets/3c8da88e-826a-43b5-85ee-6477c8494ffe" />


**Principle**: Profiles the transcriptome by sequencing cDNA derived from RNA. Captures isoforms, splicing, and expression levels.  

**Pipeline**:
1. RNA extraction & quality check (RIN score).  
2. **Library prep**: poly(A) selection or rRNA depletion; RNA fragmented and reverse-transcribed into cDNA; adapters added.  
3. Sequencing (usually Illumina).  
4. **QC & trimming** (FastQC).  
5. **Alignment** to reference genome (STAR, HISAT2) or transcriptome (Salmon, Kallisto).  
6. **Quantification** of expression (counts/TPM/FPKM).  
7. Downstream: differential expression (DESeq2, edgeR), isoform analysis, pathway enrichment.  

**Table – RNA-seq Overview**

| Feature        | Details                                                                 |
|----------------|-------------------------------------------------------------------------|
| **Principle**  | Sequence cDNA derived from RNA (polyA selected or rRNA-depleted)        |
| **Output**     | Expression profiles (counts, TPM, FPKM); isoforms; fusion transcripts   |
| **Applications** | Differential expression, splicing analysis, gene fusion discovery     |
| **Strengths**  | Genome-wide transcriptome profiling, detects novel isoforms             |
| **Limitations**| Sensitive to RNA quality, batch effects, alignment complexity           |

**RNA-Seq** = sequencing of the complete set of transcripts in a cell (the transcriptome).
-Can include:
mRNAs
non-coding RNAs
small RNAs
-Applications:
Catalogue all species of transcripts.
Determine transcriptional structure of genes (start/end sites, splicing, post-transcriptional modifications).
Quantify expression levels across conditions or developmental stages.

*Typical RNA-Seq Workflow*
Extract RNA.
Convert RNA to cDNA libraries (fragmentation + adaptor ligation).
(Optional) Amplification.
Sequence fragments with NGS (Illumina, Ion, etc.).
Reads are aligned either to:
Reference genome
Reference transcriptome
Or assembled de novo (if no reference available).
Typical read length: 30–400 bp.
Reads can be:
Exonic reads
Junction reads (span exon–exon boundaries)
Poly(A)-end reads

*Data Analysis*
Align reads to reference or assemble de novo.
Count reads overlapping each transcript → expression quantification.
Highly expressed transcripts = many mapped reads.

*Expression Metrics*
To normalize for sequencing depth and gene length:

RPKM = Reads Per Kilobase of transcript per Million mapped reads.
FPKM = Fragments Per Kilobase of transcript per Million mapped reads.
TPM = Transcripts Per Kilobase of transcript per Million mapped reads.

*Strategies for Transcript Reconstruction*
Two main approaches:
-Align-then-assemble
Reads are first aligned to reference genome.
Transcripts reconstructed from splice junctions and alignments.

-Assemble-then-align
Reads are first assembled de novo into transcripts.
Assembled transcripts are then mapped back to genome.
Works better for abundant transcripts.

*RNA-Seq for Gene Expression Analysis*
Assumption: transcript abundance is proportional to number of reads mapped.
Pipeline:
Map reads to genome.
Count reads per gene (counting matrix).
Normalize data.
Different normalization methods may lead to different results.

*Summary*
RNA-Seq enables quantification and discovery of transcripts (mRNA, ncRNA, small RNAs).
Two main strategies for transcript reconstruction: align-then-assemble and assemble-then-align.
Read counts reflect expression level → normalized using RPKM, FPKM, TPM.
Provides powerful insights into gene expression regulation, splicing, and transcriptome complexity.

---

### Genome-Wide Association Studies (GWAS)

<img width="360" height="392" alt="Screenshot 2025-08-29 160247" src="https://github.com/user-attachments/assets/bc9a7333-4c5d-4510-9ccd-2e1788e8ad01" />

**Principle**: Detects statistical associations between SNPs and phenotypes across large populations.  

**Pipeline**:
1. **Genotyping & QC** (call rate, Hardy–Weinberg equilibrium, minor allele frequency).  
2. **Population structure** correction (PCA/MDS, kinship matrices).  
3. **Association testing** (linear/logistic regression, mixed models).  
4. **Multiple testing correction** (Bonferroni, FDR).  
5. **Visualization** → Manhattan plot, QQ plot.  
6. **Interpretation** → fine-mapping, candidate gene identification, pathway analysis.  

**Manhattan plots**: each SNP plotted by genomic position (x-axis) vs –log₁₀(p) (y-axis); significant peaks highlight candidate loci.  

<img width="316" height="354" alt="Screenshot 2025-08-29 160155" src="https://github.com/user-attachments/assets/b6972539-40e2-40e7-b328-66f62025288a" />

**QQ PLOT**

<img width="353" height="333" alt="Screenshot 2025-08-29 160229" src="https://github.com/user-attachments/assets/6fb423ed-165b-4ade-a0f7-570f7a062b53" />


**Bonferroni correction**: very conservative (α/N); often supplemented by FDR to retain power.  

**Challenges**:  
- Biased toward common variants with moderate effect.  
- Rare variants require sequencing and gene-based burden tests.  
- Population stratification can inflate false positives.  
- GWAS often limited to European cohorts (ancestry bias).  

**Links to Population Genomics**:  
- **LD**: fundamental to detect associations (tag SNPs).  
- **ROH / inbreeding**: impact homozygosity and recombination patterns.
  
    <img width="319" height="334" alt="Screenshot 2025-08-29 160137" src="https://github.com/user-attachments/assets/9ae8a431-68ad-4afd-b5ec-37e095602439" />

- **PCA/MDS**: visualize and correct for population structure.  
- **Effective population size**: influences LD decay and GWAS resolution.  

**Table – GWAS Overview**

| Feature        | Details                                                                 |
|----------------|-------------------------------------------------------------------------|
| **Principle**  | Association study testing SNPs across genome vs phenotype               |
| **Pipeline**   | QC → structure correction → association → multiple testing → visualization → interpretation |
| **Output**     | Manhattan plots, QQ plots, fine-mapped loci                             |
| **Strengths**  | Genome-wide, unbiased, powerful for common variants                     |
| **Limitations**| Rare variants hard to detect, ancestry bias, population stratification issues |

---
### ROH (Runs of Homozygosity) and Inbreeding in Population Genomics

- **Runs of Homozygosity (ROH)**  
  - Long consecutive stretches of homozygous genotypes across the genome.  
  - Indicate **autozygosity**: inheritance of identical chromosomal segments from both parents.  
  - The **length and number of ROH** inform about population history:  
    - Long ROH → recent inbreeding (parents closely related).  
    - Shorter, numerous ROH → ancient bottlenecks or small effective population size.  

- **Inbreeding coefficient (F)**  
  - Probability that two alleles at a locus are **identical by descent (IBD)**.  
  - Can be estimated:  
    - From **pedigrees** (theoretical).  
    - From **genomic data** (fraction of genome in ROH, excess homozygosity).  
  - High F → more homozygosity, reduced genetic diversity, potential inbreeding depression.  

- **Population genomics perspective**  
  - ROH patterns reflect **demographic history**, **population size**, and **mating structure**.  
  - Useful for:  
    - Detecting recent **consanguinity** in individuals.  
    - Understanding **domestication** and **breeding programs** (livestock, crops).  
    - Studying endangered species with low genetic diversity.  
  - In large, outbred populations → fewer and shorter ROH.  
  - In small or isolated populations → more and longer ROH.  

**Key takeaway:**  
ROH and inbreeding coefficients are powerful indicators of both **individual-level inbreeding** and **population-level demographic history**, bridging genetics and evolutionary biology.

---
**In summary**, specialized methods (aCGH, Pool-seq, targeted panels, **WES**, methyl-seq, RNA-seq) extend genomic analysis beyond WGS. **GWAS** builds on these data to link genetic variation with traits, but requires careful handling of QC, population structure, and multiple testing to produce biologically valid results.

