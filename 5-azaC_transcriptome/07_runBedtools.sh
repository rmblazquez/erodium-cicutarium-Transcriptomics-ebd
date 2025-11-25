# Change the path to analyze different sets of BAM files

projectDir="." #"/path/to/project/directory/"
cd $projectDir/results

# Change the pattern to get different patterns of BAM files
# Left by default to generate counts for 100 kb bins
for i in $(ls star/*/*.sortedByCoord.out.bam | cut -d$'/' -f3 | cut -d$'.' -f1); do 
	bedtools coverage \
	-a $projectDir/data/genome/ecic_genome_bin100k.bed \
	-b star/$i.STAR/$i.STARAligned.sortedByCoord.out.bam -counts -sorted > bedtools/$i.bedgraph; 
done
# Use this script to generate the input for the Circos graph in Figure 1
