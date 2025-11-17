#!/usr/bin/env bash

#### SETTING THE PATH ####
# Define the path of the project
projectDir="/path/to/project/directory"
cd $projectDir

#### FILTER READS AND NUCLEOTIDES BY QUALITY WITH TRIMMOMATIC ####
# Define a path for fata file with Illumina adapters 
# Different adapters are used depending on type of sequencing machine
trimm="/path/to/file/adapters.fa" # Illumina TruSeq adapters
# Run Trimmomatic for all the fastq files
cat $projectDir/src/ecic_basenames.txt | while read line; do
  TrimmomaticPE \
  -phred33 \
  $projectDir/data/raw-seq/$line\_1.fastq.gz \
  $projectDir/data/raw-seq/$line\_2.fastq.gz \
  $projectDir/results/trimmomatic/$line\_1.trimmed.paired.fastq.gz \
  $projectDir/results/trimmomatic/$line\_1.trimmed.unpaired.fastq.gz \
  $projectDir/results/trimmomatic/$line\_2.trimmed.paired.fastq.gz \
  $projectDir/results/trimmomatic/$line\_2.trimmed.unpaired.fastq.gz \
  ILLUMINACLIP:$trimm:2:30:10 \
  LEADING:28 TRAILING:28 \
  SLIDINGWINDOW:4:15 MINLEN:30;
  mkdir $projectDir/results/fastqc-trimseq/$line\_1.trimmed.paired.fastqc;
  fastqc -o $projectDir/results/fastqc-trimseq/$line\_1.trimmed.paired.fastqc \
  $projectDir/results/trimmomatic/$line\_1.trimmed.paired.fastq.gz;
  mkdir $projectDir/results/fastqc-trimseq/$line\_2.trimmed.paired.fastqc;
  fastqc -o $projectDir/results/fastqc-trimseq/$line\_2.trimmed.paired.fastqc \
  $projectDir/results/trimmomatic/$line\_2.trimmed.paired.fastq.gz;
done

# Generate a HTML quality chek report for trimmed reads with MultiQC
cd $projectDir/results/fastqc-trimseq
multiqc $projectDir/results/fastqc-trimseq/*/*.zip
