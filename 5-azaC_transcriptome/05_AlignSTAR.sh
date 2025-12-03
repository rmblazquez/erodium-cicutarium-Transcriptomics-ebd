#!/usr/bin/env bash

#### SETTING THE PROJECT ####
# Define the path of the project
projectDir="." #"/path/to/project/directory"
cd $projectDir

#### ALIGN THE READS TO BIMP REFERENCE GENOME WITH STAR ####
cat $projectDir/src/SRR_accessions.txt | while read line; do
mkdir -p $projectDir/results/star/$line.STAR;
STAR --runThreadN 4 \
  --genomeDir $projectDir/data/genome/STAR-2.4.2a-R64_Ecic_Index \
  --readFilesIn $projectDir/results/trimmomatic/$line\_1.trimmed.paired.fastq.gz $projectDir/results/trimmomatic/$line\_2.trimmed.paired.fastq.gz \
  --sjdbGTFfile $projectDir/data/genome/GCA_046563565.1_EBD_Ecic_2.1_genomic.gtf \
  --outFileNamePrefix $projectDir/results/star/$line.STAR/$line.STAR \
  --readFilesCommand gunzip -c \
  --quantMode TranscriptomeSAM GeneCounts \
  --twopassMode Basic \
  --outSAMunmapped Within  \
  --outSAMtype BAM SortedByCoordinate;
 # Output your alignment file in BAM format, sorted by coordinate
 # Give your output files a UNIQUE name/prefix (${line} also works)
 # Gene counting for generating gene count tables
 # Two-pass mode for refining splice junctions
done

#### GENERATE MULTIQC REPORT ####
cd $projectDir/results/star/
multiqc $projectDir/results/star/
