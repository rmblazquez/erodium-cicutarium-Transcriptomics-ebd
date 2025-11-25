#!/usr/bin/env bash

#### SETTING THE PROJECT FILES ####
# Define the path of the project
projectDir="." #"/path/to/project/directory"
cd $projectDir
# Generate all subfolders
mkdir -p $projectDir/src $projectDir/data/genome $projectDir/data/raw-seq \
$projectDir/results/fastqc-rawseq $projectDir/results/fastqc-trimseq  \
$projectDir/results/htseq-count  $projectDir/results/star \
$projectDir/results/trimmomatic $projectDir/results/featureCounts \
$projectDir/results/salmon

# download genome data
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/046/563/565/GCA_046563565.1_EBD_Ecic_2.1/GCA_046563565.1_EBD_Ecic_2.1_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/046/563/565/GCA_046563565.1_EBD_Ecic_2.1/GCA_046563565.1_EBD_Ecic_2.1_genomic.gbff.gz
mv GCA_046563565.1_EBD_Ecic_2.1_genomic* $projectDir/data/genome

# download Illumina PE reads
# Generate the SRR_accessions.txt file from Table 2 from the source manuscript
cat $projectDir/src/SRR_accessions.txt | while read line; do
  fasterq-dump $line --split-files -O $projectDir/data/raw-seq/$line
done 
