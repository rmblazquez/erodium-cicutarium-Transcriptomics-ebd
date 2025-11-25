#!/usr/bin/env bash

#### SETTING THE PROJECT ####
# Define the path of the project
projectDir="." #"/path/to/project/directory"
cd $projectDir

#### GET MAPPED READ COUNTS WITH HTSEQ-COUNT (B. IMPATIENS) ####
cat $projectDir/src/SRR_accessions.txt | while read line; do
  mkdir -p $projectDir/results/htseq-count/$line.Counts;
  htseq-count -q -f bam -t exon -i gene -s no $projectDir/results/star/$line.STAR/$line.STARAligned.sortedByCoord.out.bam $projectDir/data/genome/GCA_046563565.1_EBD_Ecic_2.1_genomic.gtf > $projectDir/results/htseq-count/$line.Counts/$line.htseq-count;
done
# Count only exon overlapping counts in gene features

#### GENERATE MULTIQC REPORT ####
cd $projectDir/results/htseq-count/
multiqc $projectDir/results/htseq-count/
