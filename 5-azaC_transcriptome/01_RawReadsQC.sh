#!/usr/bin/env bash

#### SETTING THE PATH ####
# Define the path of the project
projectDir="." #"/path/to/project/directory"
cd $projectDir
# Generate a list of file names from the project
# Previously change the accession numbers to the file names in table 2 from the manuscript, make sure they end in ".fastq.gz"
ls $projectDir/data/raw-seq/*.fastq.gz | sed 's/.fastq.gz//g' $projectDir/src/ecic_basenames.txt | sort | uniq > $projectDir/src/ecic_basenames.txt

#### RUN QUALITY CHECK IN RAW SEQUENCES WITH FASTQC ####
# Run FastQC for all the project files
cat $projectDir/src/ecic_basenames.txt | while read line; do
	mkdir $projectDir/results/fastqc-rawseq/$line.fastqc
	fastqc -o $projectDir/results/fastqc-rawseq/$line.fastqc \
		$projectDir/data/raw-seq/$line.fastq.gz;
done

# Generate a HTML quality chek report for raw reads with MultiQC
cd $projectDir/results/fastq-rawseq
multiqc $projectDir/results/fastqc-rawseq/*/*.zip
