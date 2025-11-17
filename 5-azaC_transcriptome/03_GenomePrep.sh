#!/usr/bin/env bash

#### SETTING THE PATH ####
# Define the path of the project
projectDir="/path/to/project/directory/"
cd $projectDir

#### EDIT GENOME FILES TO FIT THE PIPELINE ####
gunzip $projectDir/data/genome/GCA_046563565.1_EBD_Ecic_2.1_genomic.fna.gz $projectDir/data/genome/GCA_046563565.1_EBD_Ecic_2.1_genomic.gbff.gz
grep -v '#' $projectDir/data/genome/GCA_046563565.1_EBD_Ecic_2.1_genomic.gbff > $projectDir/data/genome/GCA_046563565.1_EBD_Ecic_2.1_genomic.tmp.gbff

#### CONVERT GFF FILES INTO GTF FILES ####
# Use the Rscript gff2gtf_rtracklayer.R to convert the GFF into GTF
# NEED TO USE R version >= 3.5.0
/usr/local/bin/./R --slave --no-restore --file=$projectDir/src/gff2gtf_rtracklayer.R \
 $projectDir/data/genome/GCA_046563565.1_EBD_Ecic_2.1_genomic.tmp.gbff \
 $projectDir/data/genome/GCA_046563565.1_EBD_Ecic_2.1_genomic.tmp.gbff \
 putEGin_gene_id
