#!/usr/bin/env bash

#### SETTING THE PROJECT ####
# Define the path of the project
projectDir="." #"/path/to/project/directory"
cd $projectDir

#### ALIGN THE READS TO BIMP REFERENCE GENOME WITH STAR ####
cat $projectDir/src/basenames.txt | while read line; do
mkdir -p $projectDir/results/star/$line.STAR;
STAR --runThreadN 4 \
  --genomeDir $projectDir/data/genome/STAR-2.4.2a-R64_Ecic_Index \
  --readFilesIn $projectDir/results/trimmomatic/$line\_1.trimmed.paired.fastq.gz $projectDir/results/trimmomatic/$line\_2.trimmed.paired.fastq.gz \
  --readFilesCommand gunzip -c \
  --sjdbGTFfile $projectDir/data/genome/GCA_046563565.1_EBD_Ecic_2.1_genomic.gtf \
  --sjdbOverhang 99 \
  --limitGenomeGenerateRAM 5900000000 \
  --outSAMtype BAM SortedByCoordinate \
  --outFileNamePrefix $projectDir/results/star/$line.STAR/$line.STAR \
  --quantMode GeneCounts \
  --twopassMode Basic;
 # Do not let STAR use more than 5.9 GB of memory
 # Output your alignment file in BAM format, sorted by coordinate
 # Give your output files a UNIQUE name/prefix (${line} also works)
 # Gene counting for generating gene count tables
 # Two-pass mode for refining splice junctions
done

#### GENERATE MULTIQC REPORT ####
cd $projectDir/results/star/
multiqc $projectDir/results/star/
