# erodium-cicutarium-Transcriptomics-ebd

Author: Rubén Martín-Blázquez (@rmblazquez)

This repository contains the scripts (bash and R) and instructions to analyze RNA-seq from published data from *Erodium cicutarium*. 

## 5-azaC_transcriptome

Contains the scripts necessary to perform the analyses described in "Understanding epigenetic regulation in non-model plants: Transcriptomic responses to seed demethylation in leaves and roots of the annual herb *Erodium cicutarium* (Geraniaceae)".

Start by creating a folder called "src" in your path, and clone or copy the contents of the directory "5-azaC_transcriptome" there. The scripts contain a variable called *projectDir* with the pattern "/path/to/project/directory" that is supposed to include your path, thus they have to be modified replacing this pattern with the actual working directory's path before running them. In order to replicate the analysis, they have to be ran in such numeric order (i. e., start with 00_setDirs.sh, then 01_RawReadsQC.sh, and so on). There are auxiliary files that are not numbered, leave them be as they are, they will be required by some of the scripts.

### 00_setDirs.sh

This script will set the directory structure neede for the analyses, and will download both the draft genome and Illumina PE reads from the corresponding NCBI repositories. You will need to install **SRA toolkit** before running this script. The script uses the auxiliary file *SRR_accessions.txt*.

### 01_RawReadsQC.sh

This script will perform a quality check analysis of the raw reads. You will need to install **FastQC** and **MultiQC**before running this script.

### 02_TrimReads.sh

This script will filter low quality positions, remove adapters, and filter short sequences from the raw reads. In addition, it will run another QC with the trimmed reads. You will need to install **Trimmomatic**, **FastQC** and **MultiQC** before running this script. The script uses the auxiliary file *adapters.fa*.

### 03_GenomePrep.sh

This script will format the genome sequence and genome annotation files for the pipeline. The script uses the auxiliary file *gff2gtf_rtracklayer.R*, which requires **R version >= 3.5.0** in order to be run.
