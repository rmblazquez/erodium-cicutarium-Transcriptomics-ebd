# change the path to your working directory
setwd("/path/to/project/directory/data/genome")

library(BioCircos)

#### Genome configuration ####
# use a list with N elements (N = number of scaffolds)
# each element has an structure X = L, where
# X = scaffold name; L = scaffold length

genomeLengths <- read.table("scaffold_lengths.txt", sep = '\t', header = F)
genomeLengths <- as.data.frame(cbind(genomeLengths$V2, genomeLengths$V1))
colnames(genomeLengths) <- c("name", "length")
genomeLengths$length <- as.numeric(genomeLengths$length)
genomeLengths <- genomeLengths[order(genomeLengths$length, decreasing = T),] ###

myGenome <- list()
i <- 1
while (i <= length(genomeLengths$length)) {
  myGenome[[i]] <- as.numeric(genomeLengths$length[i])
  i <- i + 1
}
names(myGenome) <- genomeLengths$name

#### Generate 100k bins from reference genome ####
# # run the following in a Linux terminal:
# # generate BED file with 100k sized windows/bins from the reference genome
# bedtools makewindows -g ecic_genome_lengths.txt -w 100000 > ecic_genome_bin100k.bed

#### Gene expression tracks ####
# # # run the following in a Linux terminal:
# # # calculate coverage value for each 100k intervals
# # bedtools coverage -a ecic_genome_bin100k.bed -b CODE.bam -counts -sorted > coverage.bedgraph 
# # # add the expression from all treatment + tissue combinations
#### TE density track ####
# # # run in a Linux terminal
# # # compute the number of TEs in each 100k bin
# # bedtools intersect -a ecic_genome_bin100k.bed -b CicuSeq2_assembly.fasta.out.gff -c > ecic_TEs_bin100k.bed
#### Gene density track ####
# # # run in a Linux terminal
# # # compute the number of genes in each 100k bin
# # bedtools intersect -a ecic_genome_bin100k.bed -b ecic_annotation.genes.gff3 -c > ecic_genes_bin100k.bed
# # Compile all results in the file "ecic_genome_tes_genes_expression.txt"

genomeData <- read.table("ecic_genome_tes_genes_expression_v2.txt", sep = '\t', header = T)

# sort genome bins by contig length
j <- 1
k <- 1
scLengths <- c()
while (j <= length(genomeLengths$name)) {
  while (k <= length(genomeData$scaffold)) {
    if (genomeData$scaffold[k] == genomeLengths$name[j]) {
      scLengths <- c(scLengths, genomeLengths$length[j])
      k <- k + 1
      j <- 1
    } else {
      j <- j + 1
    }
  }
} # stop manually, the loop seems to be endless... ###

genomeData$lengths <- as.numeric(scLengths) ###
genomeData <- genomeData[order(genomeData$lengths, decreasing = T), ] ###

# Check for local maximums in the contig length distribution to get a cut
# N50 = 1953031
hist(genomeLengths$length[which(genomeLengths$length >= 1953031)])

# There are more 6Mb contifgs than 5 Mb contigs: use 6Mb as threshold
genomeData.bc <- genomeData[which(genomeData$lengths >= 6000000), ] ###
myGenome.bc <- myGenome[which(myGenome >= 6000000)] # reduce the number of contigs to those longer than N50

#### Configure Circos graph tracks ####
chrVert <- genomeData.bc$scaffold
chrCoor1 <- genomeData.bc$start
chrCoor2 <- genomeData.bc$end
posVert.TEs <- genomeData.bc$Tes
posVert.genes <- genomeData.bc$genes
posVert.exCL <- genomeData.bc$expr_C_L
posVert.exAL <- genomeData.bc$expr_A_L
posVert.exCR <- genomeData.bc$expr_C_R
posVert.exAR <- genomeData.bc$expr_A_R
posVert.exOC <- genomeData.bc$expr_C_C
posVert.exOA <- genomeData.bc$expr_C_A

tracks <- BioCircosBarTrack('TE density', as.character(chrVert), chrCoor1, chrCoor2, 
                             minRadius = 0.1, maxRadius = 0.34, 
                             values = posVert.TEs, color = "steelblue")
tracks <- tracks + BioCircosBackgroundTrack('TE density', minRadius = 0.1, maxRadius = 0.34,
                                            borderColors = "#FFFFFF", fillColor = "#ECECEC", borderSize = 0)
tracks <- tracks + BioCircosBarTrack('gene density', as.character(chrVert), chrCoor1, chrCoor2,
                                      minRadius = 0.35, maxRadius = 0.59,
                                      values = posVert.genes, color = "tomato")
tracks <- tracks + BioCircosBackgroundTrack('gene density', minRadius = 0.35, maxRadius = 0.59,
                                            borderColors = "#FFFFFF", fillColor = "#ECECEC", borderSize = 0)
tracks <- tracks + BioCircosBarTrack('expression CL', as.character(chrVert), chrCoor1, chrCoor2,
                                      minRadius = 0.73, maxRadius = 0.85,
                                      values = log2(posVert.exCL + 1), color = "lightgreen")
tracks <- tracks + BioCircosBackgroundTrack('expression CL', minRadius = 0.73, maxRadius = 0.85,
                                            borderColors = "#FFFFFF", fillColor = "#ECECEC", borderSize = 0)
tracks <- tracks + BioCircosBarTrack('expression CR', as.character(chrVert), chrCoor1, chrCoor2,
                                      minRadius = 0.6, maxRadius = 0.72,
                                      values = log2(posVert.exCR + 1), color = "darkorange")
tracks <- tracks + BioCircosBackgroundTrack('expression CR', minRadius = 0.6, maxRadius = 0.72,
                                            borderColors = "#FFFFFF", fillColor = "#ECECEC", borderSize = 0)
tracks <- tracks + BioCircosBarTrack('expression OC', as.character(chrVert), chrCoor1, chrCoor2,
                                      minRadius = 0.86, maxRadius = 0.99,
                                      values = log2(posVert.exOC + 1), color = "darkgreen")
tracks <- tracks + BioCircosBackgroundTrack('expression OC', minRadius = 0.86, maxRadius = 0.99,
                                            borderColors = "#FFFFFF", fillColor = "#ECECEC", borderSize = 0)

BioCircos(tracks, genome = myGenome.bc, genomeFillColor = "black",#"Spectral",
          chrPad = 0.02, displayGenomeBorder = T, 
          genomeBorderColor = "black", genomeBorderSize = 0.1,
          genomeTicksDisplay = FALSE, genomeLabelTextSize = 10,
          genomeLabelDy = 0)

#### Correlation between RE density and gene expression ####

library(Hmisc)

# Expression from all libraries
expr.all <- as.data.frame(cbind(genomeData$expr_C_L,
                                genomeData$expr_A_L,
                                genomeData$expr_C_R,
                                genomeData$expr_A_R,
                                genomeData$expr_C_C,
                                genomeData$expr_C_A))

# Correlation between gene and RE densities
rcorr(genomeData$Tes, genomeData$genes)
plot(genomeData$Tes, genomeData$genes)
lines(predict(lm(genomeData$Tes ~ genomeData$genes)), col = 'steelblue')

# Correlation between expression and RE density
rcorr(genomeData$Tes, rowSums(log(expr.all + 1))) # R = 0.47, p < 0.001, BF < 0.001
plot(genomeData$Tes, rowSums(log(expr.all + 1)))
lines(predict(lm(genomeData$Tes ~ rowSums(log(expr.all + 1)))), col = 'steelblue')

# Separated RNA-seq library: REs, Bonferroni = P-value * 8800
rcorr(genomeData$Tes, log(genomeData$expr_C_L + 1)) # R = 0.44, p < 0.001, BF < 0.001
rcorr(genomeData$Tes, log(genomeData$expr_A_L + 1)) # R = 0.44, p < 0.001, BF < 0.001
rcorr(genomeData$Tes, log(genomeData$expr_C_R + 1)) # R = 0.47, p < 0.001, BF < 0.001
rcorr(genomeData$Tes, log(genomeData$expr_A_R + 1)) # R = 0.47, p < 0.001, BF < 0.001
rcorr(genomeData$Tes, log(genomeData$expr_C_C + 1)) # R = 0.46, p < 0.001, BF < 0.001
rcorr(genomeData$Tes, log(genomeData$expr_C_A + 1)) # R = 0.46, p < 0.001, BF < 0.001

# Correlation between expression and gene density
rcorr(genomeData$genes, rowSums(log(expr.all + 1))) # R = 0.61, p < 0.001
plot(genomeData$genes, rowSums(log(expr.all + 1)))
lines(predict(lm(genomeData$genes ~ rowSums(log(expr.all + 1)))), col='tomato')

# Separated RNA-seq library: genes
rcorr(genomeData$genes, log(genomeData$expr_C_L + 1)) # R = 0.58, p < 0.001, BF < 0.001
rcorr(genomeData$genes, log(genomeData$expr_A_L + 1)) # R = 0.59, p < 0.001, BF < 0.001
rcorr(genomeData$genes, log(genomeData$expr_C_R + 1)) # R = 0.62, p < 0.001, BF < 0.001
rcorr(genomeData$genes, log(genomeData$expr_A_R + 1)) # R = 0.62, p < 0.001, BF < 0.001
rcorr(genomeData$genes, log(genomeData$expr_C_C + 1)) # R = 0.60, p < 0.001, BF < 0.001
rcorr(genomeData$genes, log(genomeData$expr_C_A + 1)) # R = 0.60, p < 0.001, BF < 0.001
