# MAIRIS: Genomic Pipeline for Haplotype-Based Multi-Gene Analysis in Diploid Organisms: Validation Using *Musca domestica* Insecticide Resistance Genes

### MAIRIS: Multipurpose Analysis of Insecticide Resistance-associated Identified SNPs, MNPs and InDels

MAIRIS is a pipeline designed for analyzing variants in protein-coding genes of diploid organisms, with a particular
focus on insecticide resistance-associated mutations. It identifies both genomic DNA variants and their resulting amino
acid changes.

## Docker Environment

### Prerequisites

- Docker
- Docker Compose
- Git

### Directory structure in Docker HOST

```
/mairis2025/         Cloned folder using git
├── houseflyVGSC/    Folder for housefly sodium channel analysis
│   └── ref          Reference files
├── houseflyAChE/    Folder for housefly acetylcholinesterase analysis
│   └── ref          Reference files
├── fastaq           Folder for copying paired-end data of NGS (including VGSC and AChE)
└── scripts          Pipeline scripts (executable in the docker container)
```

### Quick Start Docker container

1. Clone the repository in your :

```bash
git clone https://github.com/entlabo3/mairis2025.git
```

1. Build and start the Docker container:

```bash
docker-compose up --build -d
```

1. Enter the container:

```bash
docker exec -it myris2025-mairis /bin/bash
```

* Caution * The container name may be automatically modified by the program. For example, a number such as '-1' may be appended to the end. Please check the name of the running container.

Alternative version with more detail:

The container name may be automatically changed by Docker. For instance, a suffix like '-1' might be added to the container name. Please verify the actual running container name before proceeding.

You can check the running container name using:
```bash
docker ps
```


### Directory Structure in Docker container


Folders and files created or copied within the cloned folder
("/opt/mairis" in Docker container corresponds to "/mairis2025/" in Docker host)
on the Docker host are accessible from inside the container.
Therefore, you don't need to copy your data into the container.

Files required for analysis, other than those in "/opt/mairis",
are automatically configured by the Dockerfile and docker-compose.yml.

**Additional explanation for the mapping:**

```
("/opt/mairis" in Docker container corresponds to "/mairis2025/" in Docker host)
```

This configuration is established through the volume mapping in docker-compose.yml:

```yaml
volumes:
  - ./:/opt/mairis  # Mount current directory to /opt/mairis in the container
```

This means that the current directory on your host machine (where you cloned the repository) is mounted as `/opt/mairis` inside the Docker container, creating a seamless connection between your host filesystem and the container environment.

**Performance Note:**
If you have a large number of samples and sufficient Docker container resources, placing your data directly inside the container may provide faster processing speeds. However, the current volume mounting approach offers greater flexibility for data management and easier access to results.

The volume mounting approach is recommended for most users as it:
- Eliminates the need to copy large datasets into the container
- Allows real-time access to results from the host system  
- Simplifies data management and backup procedures
- Maintains data persistence even if the container is removed




```
/
└── opt
　    ├── conda         automatically installed
 　   ├── mairis        Connect to Docker host with Docker Compose volumes
 　   │    ├── fastaq     folder for your paired-end NGS data (including both VGSC and AChE)
  　  │    ├── houseflyAChE
  　  │    │    └── ref   Reference files for houseflyAChE
   　 │    ├── houseflyVGSC
   　 │    │    └── ref   Reference files for houseflyVGSC
  　  │    └── scripts       Pipeline scripts (executable from anywhere in the container)
  　  └── temp               folder for temporary files is created automatically

```

### Gene-Specific Folder Structure

Create one folder for each gene of interest (e.g., VGSC and AChE). Each folder must contain a `ref` subdirectory with gene-specific information. The `ref` folder requires a minimum of three essential files:

#### Required Reference Files

1. **Reference genome sequence**: The genomic DNA sequence of your target gene from the organism you are analyzing
2. **GFF annotation file**: A GFF3 file containing annotation information corresponding to the reference genome sequence
3. **Standard amino acid reference sequence**: An amino acid sequence that serves as the positional standard for variant annotation

##### Example Reference Files

To use as a reference for aligning short read sequence data,
we need the genomic DNA sequence and the GFF file that records the annotation information for that FASTA file.

This GFF3 file includes a custom attribute termed 'exon_number' to denote the absolute position of each CDS within the
transcript. This numbering system fulfills two functions:

1. Ordering of alternatively spliced exons: When mutually exclusive exons are present (e.g., Exon16c and Exon16d), they
   are designated with the same exon number (e.g., 16) to reflect their equivalent position within the transcript.


1. Accommodating newly discovered exons: The numbering system incorporates decimal values to facilitate future
   discoveries. For example, the discovery of a new exon between Exon09 and Exon10 can result in its assignment as
   exon_number=9.5, facilitating its integration into the current numbering system without disturbing the established
   sequence.

This custom annotation facilitates a clear comprehension of exon arrangement and alternative splicing patterns within
the structure of the target gene.

GFF3 format example:

```
##gff-version 3								
mdVGSCgenome	.	gene	1000	177517	.	+	.	ID=gene;biotype=protein_coding;Name=VGSC
mdVGSCgenome	.	transcript	1000	177517	.	+	.	ID=mRNAck;Parent=gene;Name=mRNAck;biotype=protein_coding
mdVGSCgenome	.	transcript	1000	177517	.	+	.	ID=mRNAdl;Parent=gene;Name=mRNAdl;biotype=protein_coding
mdVGSCgenome	.	transcript	1000	177517	.	+	.	ID=mRNAcl;Parent=gene;Name=mRNAcl;biotype=protein_coding
mdVGSCgenome	.	transcript	1000	177517	.	+	.	ID=mRNAdk;Parent=gene;Name=mRNAdk;biotype=protein_coding
mdVGSCgenome	.	CDS	1000	1149	.	+	0	ID=CDS01;Parent=mRNAck,mRNAdl,mRNAcl,mRNAdk;Name=Exon01;exon_number=1
mdVGSCgenome	.	CDS	27767	27922	.	+	0	ID=CDS02;Parent=mRNAck,mRNAdl,mRNAcl,mRNAdk;Name=Exon02;exon_number=2
mdVGSCgenome	.	CDS	29790	29995	.	+	0	ID=CDS03;Parent=mRNAck,mRNAdl,mRNAcl,mRNAdk;Name=Exon03;exon_number=3
mdVGSCgenome	.	CDS	52602	52730	.	+	1	ID=CDS04;Parent=mRNAck,mRNAdl,mRNAcl,mRNAdk;Name=Exon04;exon_number=4
mdVGSCgenome	.	CDS	57673	57764	.	+	1	ID=CDS05;Parent=mRNAck,mRNAdl,mRNAcl,mRNAdk;Name=Exon05;exon_number=5
mdVGSCgenome	.	CDS	57858	58134	.	+	2	ID=CDS06;Parent=mRNAck,mRNAdl,mRNAcl,mRNAdk;Name=Exon06;exon_number=6
mdVGSCgenome	.	CDS	61089	61398	.	+	1	ID=CDS07;Parent=mRNAck,mRNAdl,mRNAcl,mRNAdk;Name=Exon07;exon_number=7
mdVGSCgenome	.	CDS	62242	62547	.	+	0	ID=CDS08;Parent=mRNAck,mRNAdl,mRNAcl,mRNAdk;Name=Exon08;exon_number=8
mdVGSCgenome	.	CDS	74426	74488	.	+	2	ID=CDS09;Parent=mRNAck,mRNAdl,mRNAcl,mRNAdk;Name=Exon09;exon_number=9
mdVGSCgenome	.	CDS	79057	79237	.	+	2	ID=CDS10;Parent=mRNAck,mRNAdl,mRNAcl,mRNAdk;Name=Exon10;exon_number=10
mdVGSCgenome	.	CDS	82370	82596	.	+	2	ID=CDS11;Parent=mRNAck,mRNAdl,mRNAcl,mRNAdk;Name=Exon11;exon_number=11
mdVGSCgenome	.	CDS	82678	82777	.	+	0	ID=CDS12;Parent=mRNAck,mRNAdl,mRNAcl,mRNAdk;Name=Exon12;exon_number=12
mdVGSCgenome	.	CDS	82859	82927	.	+	0	ID=CDS13;Parent=mRNAck,mRNAdl,mRNAcl,mRNAdk;Name=Exon13;exon_number=13
mdVGSCgenome	.	CDS	96360	96613	.	+	1	ID=CDS14;Parent=mRNAck,mRNAdl,mRNAcl,mRNAdk;Name=Exon14;exon_number=14
mdVGSCgenome	.	CDS	96678	96851	.	+	0	ID=CDS15;Parent=mRNAck,mRNAdl,mRNAcl,mRNAdk;Name=Exon15;exon_number=15
mdVGSCgenome	.	CDS	97063	97225	.	+	0	ID=CDS16c;Parent=mRNAck,mRNAcl;Name=Exon16c;exon_number=16
mdVGSCgenome	.	CDS	97964	98126	.	+	0	ID=CDS16d;Parent=mRNAdl,mRNAdk;Name=Exon16d;exon_number=16
mdVGSCgenome	.	CDS	103857	104044	.	+	1	ID=CDS17;Parent=mRNAck,mRNAdl,mRNAcl,mRNAdk;Name=Exon17;exon_number=17
mdVGSCgenome	.	CDS	104175	104379	.	+	2	ID=CDS18;Parent=mRNAck,mRNAdl,mRNAcl,mRNAdk;Name=Exon18;exon_number=18
mdVGSCgenome	.	CDS	114339	114520	.	+	2	ID=CDS19;Parent=mRNAck,mRNAdl,mRNAcl,mRNAdk;Name=Exon19;exon_number=19
mdVGSCgenome	.	CDS	116815	117300	.	+	0	ID=CDS20;Parent=mRNAck,mRNAdl,mRNAcl,mRNAdk;Name=Exon20;exon_number=20
mdVGSCgenome	.	CDS	118182	118355	.	+	0	ID=CDS21;Parent=mRNAck,mRNAdl,mRNAcl,mRNAdk;Name=Exon21;exon_number=21
mdVGSCgenome	.	CDS	122535	122657	.	+	0	ID=CDS22k;Parent=mRNAck,mRNAdk;Name=Exon22k;exon_number=22
mdVGSCgenome	.	CDS	132803	132925	.	+	0	ID=CDS22l;Parent=mRNAdl,mRNAcl;Name=Exon22l;exon_number=22
mdVGSCgenome	.	CDS	140158	140280	.	+	0	ID=CDS23;Parent=mRNAck,mRNAdl,mRNAcl,mRNAdk;Name=Exon23;exon_number=23
mdVGSCgenome	.	CDS	146330	146524	.	+	2	ID=CDS24;Parent=mRNAck,mRNAdl,mRNAcl,mRNAdk;Name=Exon24;exon_number=24
mdVGSCgenome	.	CDS	146896	147141	.	+	0	ID=CDS25;Parent=mRNAck,mRNAdl,mRNAcl,mRNAdk;Name=Exon25;exon_number=25
mdVGSCgenome	.	CDS	152552	152822	.	+	2	ID=CDS26;Parent=mRNAck,mRNAdl,mRNAcl,mRNAdk;Name=Exon26;exon_number=26
mdVGSCgenome	.	CDS	168530	168834	.	+	0	ID=CDS27;Parent=mRNAck,mRNAdl,mRNAcl,mRNAdk;Name=Exon27;exon_number=27
mdVGSCgenome	.	CDS	176570	177517	.	+	0	ID=CDS28;Parent=mRNAck,mRNAdl,mRNAcl,mRNAdk;Name=Exon28;exon_number=28```
```


#### Importance of Amino Acid Reference Standards

The amino acid reference sequence is particularly crucial for genes with high variability. For example, insect sodium channels exhibit extensive amino acid insertions and deletions even within the same species, making direct amino acid numbering comparisons nearly impossible. In such cases, researchers typically use the amino acid sequence from the earliest foundational papers (such as the original housefly sequence) as the standard reference for consistent numbering across studies.

#### Choosing an Appropriate Reference

For genes without established numbering conventions:
- Use sequences from major genome projects as references
- Select amino acid sequences from well-characterized standard strains or individuals that you designate as your reference
- Ensure consistency with sequences commonly used in your research field

This standardized approach ensures that variant positions can be consistently compared across different studies and populations.



## Pipeline Features

- Processes paired-end NGS data
- Handles multiple samples simultaneously
- Identifies SNPs, MNPs, insertions, and deletions
- Translates DNA variants to amino acid changes
- Performs clustering analysis of samples based on variant patterns
- Supports any diploid organism
- Applicable to any protein-coding gene sequence

## Docker Container Basics

### Base Image

[condaforge/miniforge3 (release 24.9.2):](https://hub.docker.com/r/condaforge/miniforge3) This image provides a minimal Anaconda installation with the conda package manager, ideal for creating reproducible scientific environments.

### System Packages (installed by dockerfile)

default-jre: Java Runtime Environment, likely needed for Picard. Version: 2:1.11-72

### Conda Packages  (installed by dockerfile)

* python: Programming language (v3.10).
* pandas: Data analysis library (v2.2.3).
* psutil: System monitoring library (v6.1.0).
* scikit-learn: Machine learning library (v1.5.1).
* matplotlib: Data visualization library (v3.9.1).

### Bioinformatics Tools  (installed by dockerfile)

- bbmap: Short read aligner (v37.62).
- beagle: Genotype phasing and imputation tool (v4.1).
- biopython: Library for biological computation (v1.84).
- bcftools: Tools for manipulating VCF and BCF files (v1.20).
- freebayes: Variant caller (v1.3.8).
- gffread: GFF/GTF file manipulation tool (v0.12.7).
- mafft: Multiple sequence alignment program (v7.525).
- picard: Set of tools for working with high-throughput sequencing data (v3.2.0). Note that a specific version (3.3.0) is also downloaded from GitHub.
- pyfaidx: FASTA indexing and sequence retrieval library (v0.8.1.1).
- samtools: Tools for manipulating SAM/BAM files (v1.19).
- scipy: Scientific computing library (v1.13.1).
- strobealign: Read aligner (v0.14.0).
- tabix: Tool for indexing and retrieving data from TAB-delimited files (v1.11).

## Input Files

The locations of the fastaq and ref folders can be easily customized by modifying the program configuration.



1. Paired-end FASTQ files (place in `./fastaq/`)
    - Format: `11111111111_*_R1.fastq.gz` and `11111111111*_R2.fastq.gz`
    - Our system uses the first 11 characters of fastaq filenames as unique sample identification codes. Additionally, we rely on Illumina's default filename convention with '_R1.fastq.gz' and '_R2.fastq.gz' suffixes to identify paired-end sequencing files. Code modification will be required if your files do not follow this naming convention.

2. Reference files (place in `./ref/`)
    - Reference DNA sequence for NGS results: `ref.fasta`
    - GFF3 annotation for ref: `ref.gff3`
    - Reference Aminoacid sequence:`RefAA.fasta`

## Running the Pipeline

1. Place input files in appropriate directories
2. Enter the Docker container
3. Navigate to your analysis directory
4. Run the pipeline. The command is "pipeline.sh"    (Do not use "./").

```bash
cd ./houseflyVGSC
pipeline.sh
```

## Output

The pipeline generates multiple output directories containing:

- Quality-filtered reads
- Alignment files with sorted (BAM)
- Variant calls (VCF)
- DNA and amino acid variant tables
- Hamming distance among haplotypes
- Clustering analysis results


## Pipeline Steps

1. Reference file processing
2. Sample list creation
3. Quality control
4. Read alignment
5. Duplicate marking
6. Variant calling
7. Variant phasing
8. Consensus sequence generation
9. Variant analysis
10. Hamming distance-based gene network construction

## Mapping results in JBrowse2

If you want to check the mapping results, you can use [J-Browse2](https://jbrowse.org/jb2/) the following procedure:

1. Install JBrowse2 on your PC and start it.
2. Select `Open sequence File(s)` and Input the information:
   - Assembly name : any but required
   - Type : `IndexedFastaAdapter`
   - File: `ref_trimmed.fasta` in ref folder
   - faiLocation: `ref_trimmed.fasta.fai` in ref folder
3. Select `SUBIT`, `Linear genome view`, `LAUNCH VIEW` and `SHOW ALL REAGION ASSMBLY`.
4. `OPEN TRACK SLECTOR` and check first track `Reference sequence(VGSC)`etc.
5. Select `Add track` (`+` mark at Right bottom corner or three short lines at Left upper corner ) and input the information for Exon information:
   - Main file: `ref.gff3` in ref folder
   - `NEXT` and `ADD`
6. Select `Add track` and input the information to display the mapping results:
   - Main file - File: `sampleID.sorted.markedup.bam` in the `0130_picard_markduplicates` folder
   - Index file - File: `sampleID.sorted.markedup.bam.bai` in the `0130_picard_markduplicates` folder
   - `NEXT` and `ADD`

![JBrowse2 example](fig/jbmapping.png)

## Notes

- Designed for protein-coding sequences
- Suitable for analyzing variants in diploid organisms
- Compatible with SNPs, MNPs, Insertion and Deletion variants
- Can handle multiple transcript variants
- Supports detection of both synonymous and non-synonymous mutations

## License

This software is released under the GNU General Public License v3.0.

## Citation

If you use MAIRIS in your research, please cite:
[Citation information will be added]
