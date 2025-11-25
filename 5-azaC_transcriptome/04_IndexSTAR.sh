#!/usr/bin/env bash

#### SETTING THE PROJECT ####
# Define the path of the project
projectDir="." #"/path/to/project/directory"
cd $projectDir

#### INDEX BIMP GENOME WITH STAR ####
# Create a folder for the index files
mkdir -p $projectDir/data/genome/STAR-2.4.2a-R64_Ecic_Index
# run STAR in index generation mode
STAR --runThreadN 4 \
 --runMode genomeGenerate \
 --genomeDir $projectDir/data/genome/STAR-2.4.2a-R64_Ecic_Index \
 --genomeFastaFiles $projectDir/data/genome/GCA_046563565.1_EBD_Ecic_2.1_genomic.fna.gz \
 --sjdbGTFfile $projectDir/data/genome/GCA_046563565.1_EBD_Ecic_2.1_genomic.gbff.gz \
 --sjdbOverhang 99
