# change the path to your working directory
setwd("/path/to/project/directory/results/htseq-count")

#### E. CICUTARIUM AZACITIDINE: TISSUE + AGE DATASET ####

library(DESeq2)

# Filter genes with >= 100 reads in all samples by filtering with apply(x, sum) and then selecting resulting genes
counts <- read.csv("counts-aza-age-greater100.csv", sep = ";", header = T, row.names = 1)

library(stringr)
samples <- colnames(counts)
azacytidine <- as.factor(str_sub(samples, -5, -5))
tissue <- as.factor(str_sub(samples, -3, -3))
age <- as.factor(str_sub(samples, -1,-1))
group <- as.factor(paste0(azacytidine, tissue, age))
group <- relevel(group, "CLY")

table <- data.frame(sampleName = samples, fileName = samples, condition = group)
expDesign <- model.matrix(~ 0 + group) # Define contrasts with group variable
ddsMatrix <- DESeqDataSetFromMatrix(countData = counts, colData = table, design = expDesign)

dds <- ddsMatrix[rowSums(counts(ddsMatrix)) > 1, ]
dds <- DESeq(ddsMatrix)
normCounts <- counts(dds, normalized = TRUE)
write.table(normCounts, file = "ecic-all_normcounts_deseq2.txt", 
            sep = "\t", row.names = rownames(normCounts), col.names = colnames(normCounts), 
            quote = F, append = F)

## Compare raw vs normalized counts ##

# Color vector: white = CLY, red = AYL, orange = CRY, brown = ARY, green = CLO, darkgreen = ALO
sampleColor <- c("white", "orange", "red", "brown", 
                 "white", "orange", "red", "brown", 
                 "white", "red", "brown", "white", 
                 "orange", "red", "brown", "white", 
                 "orange", "red", "brown", "white", 
                 "orange", "red", "brown", "green", 
                 "darkgreen", "green", "darkgreen", "green", 
                 "darkgreen", "darkgreen", "green", "darkgreen", 
                 "green", "darkgreen", "green")

par(mfcol = c(2,1))
boxplot(log(counts,2), col = sampleColor)
boxplot(log(normCounts, 2), col = sampleColor)

## NMDS OF TOP 500 MOST EXPRESSED GENES ##

library(edgeR)
plotMDS(dds, gene.selection = "common", top = 500, dim.plot = c(1,2), labels = NULL, 
        pch = 21,
        bg = sampleColor)

## DIFFERENTIAL EXPRESSION ANALYSIS: GENERAL EFFECTS ##

# contrast vector: group positions are groupCLY groupALO groupALY groupARY groupCLO groupCRY

# Control vs Azacytidine (all samples)
res.CvsA <- results(dds, contrast = c(1,-1,-1,-1,1,1)) # Activate ALO+ALY+ARY vs CLY+CLO+CRY comparison
data.CvsA <- data.frame(res.CvsA)
data.CvsA <- data.CvsA[with(data.CvsA, order(padj)), ]
write.table(res.CvsA, file = "18042024_deseq2_results_ecic-aza_CvsA.txt",
            sep = "\t", row.names = rownames(res.CvsA), col.names = colnames(res.CvsA),
            quote = F, append = F)
dim(data.CvsA[which(data.CvsA$padj <= 0.05),]) # 164 DEGs
dim(data.CvsA[which(data.CvsA$padj <= 0.05 & data.CvsA$log2FoldChange < 0),]) # 147 Up-regulated by azacytidine
dim(data.CvsA[which(data.CvsA$padj <= 0.05 & data.CvsA$log2FoldChange > 0),]) # 17 Down-regulated by azacytidine

# Leaf vs Root (excluded adult tissues)
res.LvsR <- results(dds, contrast = c(1,0,1,-1,0,-1)) # Activate CLY+ALY vs CRY+ARY comparison, ignore CLO+ALO
data.LvsR <- data.frame(res.LvsR)
data.LvsR <- data.LvsR[with(data.LvsR, order(padj)), ]
write.table(res.LvsR, file = "18042024_deseq2_results_ecic-tissue_LvsR.txt",
            sep = "\t", row.names = rownames(res.LvsR), col.names = colnames(res.LvsR),
            quote = F, append = F)
dim(data.LvsR[which(data.LvsR$padj <= 0.05),]) # 54424 DEGs
dim(data.LvsR[which(data.LvsR$padj <= 0.05 & data.LvsR$log2FoldChange < 0),]) # 30053 Up-regulated in root
dim(data.LvsR[which(data.LvsR$padj <= 0.05 & data.LvsR$log2FoldChange > 0),]) # 24371 Up-regulated in leaf

# Young vs Adult (excluded root samples)
res.YvsO <- results(dds, contrast = c(1,-1,1,0,-1,0)) # Activate CLY+ALY vs CLO+ALO comparison, ignore CRY+ARY
data.YvsO <- data.frame(res.YvsO)
data.YvsO <- data.YvsO[with(data.YvsO, order(padj)), ]
write.table(res.YvsO, file = "18042024_deseq2_results_ecic-age-YvsO.txt", 
            sep = "\t", row.names = rownames(res.YvsO), col.names = colnames(res.YvsO), 
            quote = F, append = F)
dim(data.YvsO[which(data.YvsO$padj <= 0.05),]) # 37102 DEGs
dim(data.YvsO[which(data.YvsO$padj <= 0.05 & data.YvsO$log2FoldChange < 0),]) # 17391 Up-regulated in Adult
dim(data.YvsO[which(data.YvsO$padj <= 0.05 & data.YvsO$log2FoldChange > 0),]) # 19711 Up-regulated in Young

# Check expression profiles from top DEG visually

DEGs.CvsA <- data.CvsA[which(data.CvsA$padj <= 0.05 & data.CvsA$log2FoldChange < 0),]
topDEG.CvsA <- rownames(DEGs.CvsA[3,])
topDEG.CvsAnormCounts <- normCounts[which(rownames(normCounts) == topDEG.CvsA),]
barplot(topDEG.CvsAnormCounts, col = sampleColor)

DEGs.LvsR <- data.LvsR[which(data.LvsR$padj <= 0.05 & data.LvsR$log2FoldChange < 0),]
topDEG.LvsR <- rownames(DEGs.LvsR[1,])
topDEG.LvsRnormCounts <- normCounts[which(rownames(normCounts) == topDEG.LvsR),]
barplot(topDEG.LvsRnormCounts, col = sampleColor)

DEGs.YvsO <- data.YvsO[which(data.YvsO$padj <= 0.05 & data.YvsO$log2FoldChange < 0),]
topDEG.YvsO <- rownames(DEGs.YvsO[1,])
topDEG.YvsOnormCounts <- normCounts[which(rownames(normCounts) == topDEG.YvsO),]
barplot(topDEG.YvsOnormCounts, col = sampleColor)

## TEST LOG FOLD CHANGE STANDARD ERROR DISTRIBUTION FROM DEGS ##

DEGs.CvsA <- data.CvsA[which(data.CvsA$padj <= 0.05),]
DEGs.CvsA$group <- c(rep("CvsA", length(DEGs.CvsA$baseMean)))
DEGs.LvsR <- data.LvsR[which(data.LvsR$padj <= 0.05),]
DEGs.LvsR$group <- c(rep("LvsR", length(DEGs.LvsR$baseMean)))
DEGs.YvsO <- data.YvsO[which(data.YvsO$padj <= 0.05),]
DEGs.YvsO$group <- c(rep("YvsO", length(DEGs.YvsO$baseMean)))

DEGs.lfcSE <- as.data.frame(c(DEGs.CvsA$lfcSE, DEGs.LvsR$lfcSE, DEGs.YvsO$lfcSE))
colnames(DEGs.lfcSE) <- "lfcSE"
DEGs.lfcSE$group <- as.factor(c(DEGs.CvsA$group, DEGs.LvsR$group, DEGs.YvsO$group))

library(ggplot2)
library(ggpubr)
library(RColorBrewer)

ggplot(DEGs.lfcSE, aes(x = lfcSE, fill = group)) +
  geom_density(alpha = 0.7) +
  labs(title = "Kernel Density Plot of DEG logFC standard errors",
       x = "SE(Logarithm Fold Change)",
       y = "Density") +
  theme_light()

ggboxplot(DEGs.lfcSE, x = "group", y = "lfcSE", color = "group", 
          repel = T, font.label = list(size = 14, face = "plain"), 
          add = "jitter", shape = "group",
          ggtheme = theme_light()) + theme(legend.position="none") +
          labs(title = "DEG logFC standard errors",
               x = "Group",
               y = "SE(Logarithm Fold Change)") 

kruskal.test(DEGs.lfcSE$lfcSE, DEGs.lfcSE$group) # K-W X2 = 1979.8, df = 2, P < 0.001
# Average + SD lfcSE: CvsA = 3.76 + 3.17, LvsR = 0.69 + 0.61. YvsO = 0.54 + 0.51

# boxplot(DEGs.CvsA$lfcSE, DEGs.LvsR$lfcSE, DEGs.YvsO$lfcSE)



## DIFFERENTIAL EXPRESSION ANALYSIS: WITHIN CATEGORY CONTRASTS ##

# contrast vector: group positions are groupCLY groupALO groupALY groupARY groupCLO groupCRY

# Azacytidine within leaf (youngs)
res.leaf <- results(dds, contrast = c(1,0,-1,0,0,0)) # Activate CLY vs ALY comparison from group variable
data.leaf <- data.frame(res.leaf)
data.leaf <- data.leaf[with(data.leaf, order(padj)), ]
write.table(res.leaf, file = "deseq2_results_ecic-aza_CLYvsALY.txt",
            sep = "\t", row.names = rownames(res.leaf), col.names = colnames(res.leaf),
            quote = F, append = F)
dim(data.leaf[which(data.leaf$padj <= 0.05),]) # 99 DEGs
dim(data.leaf[which(data.leaf$padj <= 0.05 & data.leaf$log2FoldChange < 0),]) # 59 Up-regulated by azacytidine
dim(data.leaf[which(data.leaf$padj <= 0.05 & data.leaf$log2FoldChange > 0),]) # 40 Down-regulated by azacytidine

# Azacytidine within root
res.root <- results(dds, contrast = c(0,0,0,-1,0,1)) # Activate CRY vs ARY comparison from group variable
data.root <- data.frame(res.root)
data.root <- data.root[with(data.root, order(padj)), ]
write.table(res.root, file = "deseq2_results_ecic-aza_CRYvsARY.txt",
            sep = "\t", row.names = rownames(res.root), col.names = colnames(res.root),
            quote = F, append = F)
dim(data.root[which(data.root$padj <= 0.05),]) # 487 DEGs
dim(data.root[which(data.root$padj <= 0.05 & data.root$log2FoldChange < 0),]) # 464 Up-regulated by azacytidine
dim(data.root[which(data.root$padj <= 0.05 & data.root$log2FoldChange > 0),]) # 23 Down-regulated by azacytidine

# Control leaf vs root
res.control <- results(dds, contrast = c(1,0,0,0,0,-1)) # Activate CLY vs CRY comparison from group variable
data.control <- data.frame(res.control)
data.control <- data.control[with(data.control, order(padj)), ]
write.table(res.control, file = "deseq2_results_ecic-tissue_CLYvsCRY.txt",
            sep = "\t", row.names = rownames(res.control), col.names = colnames(res.control),
            quote = F, append = F)
dim(data.control[which(data.control$padj <= 0.05),]) # 48564 DEGs
dim(data.control[which(data.control$padj <= 0.05 & data.control$log2FoldChange < 0),]) # 26812 Up-regulated in root
dim(data.control[which(data.control$padj <= 0.05 & data.control$log2FoldChange > 0),]) # 21752 Up-regulated in leaf

# Azacytidine leaf vs root
res.aza <- results(dds, contrast = c(0,0,1,-1,0,0)) # Activate ALY vs ARY comparison from group variable
data.aza <- data.frame(res.aza)
data.aza <- data.aza[with(data.aza, order(padj)), ]
write.table(res.aza, file = "deseq2_results_ecic-tissue_ALYvsARY.txt",
            sep = "\t", row.names = rownames(res.aza), col.names = colnames(res.aza),
            quote = F, append = F)
dim(data.aza[which(data.aza$padj <= 0.05),]) # 48670 DEGs
dim(data.aza[which(data.aza$padj <= 0.05 & data.aza$log2FoldChange < 0),]) # 27453 Up-regulated in root
dim(data.aza[which(data.aza$padj <= 0.05 & data.aza$log2FoldChange > 0),]) # 21217 Up-regulated in leaf

# Azacytidine within adults
res.COvsAO <- results(dds, contrast = c(0,-1,0,0,1,0)) # Activate CLO vs ALO comparison from group variable
data.COvsAO <- data.frame(res.COvsAO)
data.COvsAO <- data.COvsAO[with(data.COvsAO, order(padj)), ]
write.table(res.COvsAO, file = "deseq2_results_ecic-age-CLOvsALO.txt", 
            sep = "\t", row.names = rownames(res.COvsAO), col.names = colnames(res.COvsAO), 
            quote = F, append = F)
dim(data.COvsAO[which(data.COvsAO$padj <= 0.05),]) # 187 DEGs
dim(data.COvsAO[which(data.COvsAO$padj <= 0.05 & data.COvsAO$log2FoldChange < 0),]) # 142 Up-regulated by Azacytidine
dim(data.COvsAO[which(data.COvsAO$padj <= 0.05 & data.COvsAO$log2FoldChange > 0),]) # 45 Down-regulated by Azacytidine

# Young vs Adult (control leaf samples)
res.CYvsCO <- results(dds, contrast = c(1,0,0,0,-1,0)) # Activate CLY vs CLO comparison from group variable
data.CYvsCO <- data.frame(res.CYvsCO)
data.CYvsCO <- data.CYvsCO[with(data.CYvsCO, order(padj)), ]
write.table(res.CYvsCO, file = "deseq2_results_ecic-age-CLYvsCLO.txt", 
            sep = "\t", row.names = rownames(res.CYvsCO), col.names = colnames(res.CYvsCO), 
            quote = F, append = F)
dim(data.CYvsCO[which(data.CYvsCO$padj <= 0.05),]) # 26686 DEGs
dim(data.CYvsCO[which(data.CYvsCO$padj <= 0.05 & data.CYvsCO$log2FoldChange < 0),]) # 12635 Up-regulated in Adult
dim(data.CYvsCO[which(data.CYvsCO$padj <= 0.05 & data.CYvsCO$log2FoldChange > 0),]) # 14051 Up-regulated in Young

# Young vs Adult (azacytidine leaf samples)
res.AYvsAO <- results(dds, contrast = c(0,1,-1,0,0,0)) # Activate ALY vs ALO comparison from group variable
data.AYvsAO <- data.frame(res.AYvsAO)
data.AYvsAO <- data.AYvsAO[with(data.AYvsAO, order(padj)), ]
write.table(res.AYvsAO, file = "deseq2_results_ecic-age-ALYvsALO.txt", 
            sep = "\t", row.names = rownames(res.AYvsAO), col.names = colnames(res.AYvsAO), 
            quote = F, append = F)
dim(data.AYvsAO[which(data.AYvsAO$padj <= 0.05),]) # 28846 DEGs
dim(data.AYvsAO[which(data.AYvsAO$padj <= 0.05 & data.AYvsAO$log2FoldChange < 0),]) # 15642 Up-regulated in Adult
dim(data.AYvsAO[which(data.AYvsAO$padj <= 0.05 & data.AYvsAO$log2FoldChange > 0),]) # 13204 Up-regulated in Young


#### PLOTS ####

### Number of Up-regulated and Down-regulated DEGs with ggplot2 ###

## General effects ## 

# Azacytidine
geneExpCvsA <- cbind(c(147, 17))
colnames(geneExpCvsA) <- "Azacytidine"
rownames(geneExpCvsA) <- c("Up-regulated", "Down-regulated")

library(reshape2)
geneExpCvsA.melt <- melt(geneExpCvsA)

library(ggplot2)
plotCvsA.2 <- ggplot(data = geneExpCvsA.melt, aes(x = Var1, y = value, fill = Var1)) +  
  geom_bar(stat = "identity", position = position_dodge())  + theme_classic() +
  xlab("Treatment") + ylab("# DEGs") + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) + theme(legend.position = "none")

# Tissue
geneExpLvsR <- cbind(c(24371, 30053))
colnames(geneExpLvsR) <- "Tissue"
rownames(geneExpLvsR) <- c("Leaf", "Root")

library(reshape2)
geneExpLvsR.melt <- melt(geneExpLvsR)

library(ggplot2)
plotLvsR.2 <- ggplot(data = geneExpLvsR.melt, aes(x = Var1, y = value, fill = Var1)) +  
  geom_bar(stat = "identity", position = position_dodge())  + theme_classic() +
  xlab("Tissue") + ylab("# DEGs") + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) + theme(legend.position = "none")

# Age
geneExpYvsO <- cbind(c(19711, 17391))
colnames(geneExpYvsO) <- "Age"
rownames(geneExpYvsO) <- c("Young", "Adult")

library(reshape2)
geneExpYvsO.melt <- melt(geneExpYvsO)

library(ggplot2)
plotYvsO.2 <- ggplot(data = geneExpYvsO.melt, aes(x = Var1, y = value, fill = Var1)) +  
  geom_bar(stat = "identity", position = position_dodge())  + theme_classic() +
  xlab("Age") + ylab("# DEGs") + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) + theme(legend.position = "none")


## effects of treatments within categories

# azacytidine within Leaf (also azacytidine within young)
geneExpLCA <- cbind(c(59, 40))
colnames(geneExpLCA) <- "AzacytidineLeaf"
rownames(geneExpLCA) <- c("Up-regulated", "Down-regulated")

library(reshape2)
geneExpLCA.melt <- melt(geneExpLCA)

library(ggplot2)
plotLCA.2 <- ggplot(data = geneExpLCA.melt, aes(x = Var1, y = value, fill = Var1)) +  
  geom_bar(stat = "identity", position = position_dodge())  + theme_classic() +
  xlab("Treatment (leaf)") + ylab("# DEGs") + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) + theme(legend.position = "none")

# azacytidine within Root
geneExpRCA <- cbind(c(464, 23))
colnames(geneExpRCA) <- "AzacytidineRoot"
rownames(geneExpRCA) <- c("Up-regulated", "Down-regulated")

library(reshape2)
geneExpRCA.melt <- melt(geneExpRCA)

library(ggplot2)
plotRCA.2 <- ggplot(data = geneExpRCA.melt, aes(x = Var1, y = value, fill = Var1)) +  
  geom_bar(stat = "identity", position = position_dodge())  + theme_classic() +
  xlab("Treatment (root)") + ylab("# DEGs") + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) + theme(legend.position = "none")

# control leaf vs root
geneExpCLR <- cbind(c(21752, 26812))
colnames(geneExpCLR) <- "Control tissues"
rownames(geneExpCLR) <- c("Leaf", "Root")

library(reshape2)
geneExpCLR.melt <- melt(geneExpCLR)

library(ggplot2)
plotCLR.2 <- ggplot(data = geneExpCLR.melt, aes(x = Var1, y = value, fill = Var1)) +  
  geom_bar(stat = "identity", position = position_dodge())  + theme_classic() +
  xlab("Tissue (control)") + ylab("# DEGs") + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) + theme(legend.position = "none")

# Azacytidine leaf vs root
geneExpALR <- cbind(c(21217, 27453))
colnames(geneExpALR) <- "Azacytidine Tissues"
rownames(geneExpALR) <- c("Leaf", "Root")

library(reshape2)
geneExpALR.melt <- melt(geneExpALR)

library(ggplot2)
plotALR.2 <- ggplot(data = geneExpALR.melt, aes(x = Var1, y = value, fill = Var1)) +  
  geom_bar(stat = "identity", position = position_dodge())  + theme_classic() +
  xlab("Tissue (azacytidine)") + ylab("# DEGs") + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) + theme(legend.position = "none")

# Azacytidine within young (also azacytidine within leaf)
plotCYvsAY.2 <- ggplot(data = geneExpLCA.melt, aes(x = Var1, y = value, fill = Var1)) +  
  geom_bar(stat = "identity", position = position_dodge())  + theme_classic() +
  xlab("Treatment (young)") + ylab("# DEGs") + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) + theme(legend.position = "none")

# Azacytidine within adult
geneExpCOvsAO <- cbind(c(142, 45))
colnames(geneExpCOvsAO) <- "Azacytidine"
rownames(geneExpCOvsAO) <- c("Up-regulated", "Down-regulated")

library(reshape2)
geneExpCOvsAO.melt <- melt(geneExpCOvsAO)

library(ggplot2)
plotCOvsAO.2 <- ggplot(data = geneExpCOvsAO.melt, aes(x = Var1, y = value, fill = Var1)) +  
  geom_bar(stat = "identity", position = position_dodge())  + theme_classic() +
  xlab("Azacytidine (adult)") + ylab("# DEGs") + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) + theme(legend.position = "none")

# Age within controls
geneExpCYvsCO <- cbind(c(14051, 12635))
colnames(geneExpCYvsCO) <- "Age"
rownames(geneExpCYvsCO) <- c("Young", "Adult")

library(reshape2)
geneExpCYvsCO.melt <- melt(geneExpCYvsCO)

library(ggplot2)
plotCYvsCO.2 <- ggplot(data = geneExpCYvsCO.melt, aes(x = Var1, y = value, fill = Var1)) +  
  geom_bar(stat = "identity", position = position_dodge())  + theme_classic() +
  xlab("Age (control)") + ylab("# DEGs") + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) + theme(legend.position = "none")

# Age within azacytidine
geneExpAYvsAO <- cbind(c(13204, 15642))
colnames(geneExpAYvsAO) <- "Age"
rownames(geneExpAYvsAO) <- c("Young", "Adult")

library(reshape2)
geneExpAYvsAO.melt <- melt(geneExpAYvsAO)

library(ggplot2)
plotAYvsAO.2 <- ggplot(data = geneExpAYvsAO.melt, aes(x = Var1, y = value, fill = Var1)) +  
  geom_bar(stat = "identity", position = position_dodge())  + theme_classic() +
  xlab("Age (azacytidine)") + ylab("# DEGs") + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) + theme(legend.position = "none")

# all general effects
library(gridExtra)
grid.arrange(plotCvsA.2, plotLvsR.2, plotYvsO.2, ncol = 3, nrow = 1)
# tissue differences
grid.arrange(plotLCA.2, plotRCA.2, plotCLR.2, plotALR.2, ncol = 2, nrow = 2)
# age differences
grid.arrange(plotCYvsAY.2, plotCOvsAO.2, plotCYvsCO.2, plotAYvsAO.2, ncol = 2, nrow = 2)
# within group differences
plotRCA.3 <- plotRCA.2 + ylim(0,500) + xlab("Young root")
plotCYvsAY.3 <- plotCYvsAY.2 + ylim(0,500) + xlab("Young leaf")
plotCOvsAO.3 <- plotCOvsAO.2 + ylim(0,500) + xlab("Adult leaf")
grid.arrange(plotRCA.3, plotCYvsAY.3, plotCOvsAO.3, ncol = 3, nrow = 1)

### DEG PCAs ###

# PCA for CvsA DEGs
DEGs.CvsA <- rownames(data.CvsA[which(data.CvsA$padj <= 0.05),])
DEGs.normCounts.CvsA <- normCounts[which(DEGs.CvsA %in% rownames(normCounts)), ]

# PCA for LvsR DEGs
DEGs.LvsR <- rownames(data.LvsR[which(data.LvsR$padj <= 0.05),])
DEGs.normCounts.LvsR <- normCounts[which(DEGs.LvsR %in% rownames(normCounts)), ]

# PCA for YvsO DEGs
DEGs.YvsO <- rownames(data.YvsO[which(data.YvsO$padj <= 0.05),])
DEGs.normCounts.YvsO <- normCounts[which(DEGs.YvsO %in% rownames(normCounts)), ]

par(mfcol = c(1,3))
plotMDS(DEGs.normCounts.CvsA, gene.selection = "common", top = 500, dim.plot = c(1,2), labels = NULL, 
        pch = 21,
        bg = sampleColor)
plotMDS(DEGs.normCounts.LvsR, gene.selection = "common", top = 500, dim.plot = c(1,2), labels = NULL, 
        pch = 21,
        bg = sampleColor)
plotMDS(DEGs.normCounts.YvsO, gene.selection = "common", top = 500, dim.plot = c(1,2), labels = NULL, 
        pch = 21,
        bg = sampleColor)

### Volcano Plots ###

library(ggplot2)
library("RColorBrewer")

volcDataA <- read.table('deseq2_results_ecic-aza_CvsA.txt', sep = '\t', header = TRUE)
volcDataA$diffexpr <- "NO"
volcDataA$diffexpr[volcDataA$log2FoldChange > 0 & volcDataA$padj < 0.05] <- "DOWN"
volcDataA$diffexpr[volcDataA$log2FoldChange < 0 & volcDataA$padj < 0.05] <- "UP"
volcPlotA <- ggplot(data = volcDataA, aes(x = log2FoldChange, y = -log10(padj), colour = diffexpr)) +	
  geom_point(alpha = 0.4, size = 1.75) + 
  theme_minimal() + labs(title = "Control vs Azacytidine") +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5), axis.line = element_line(colour = "black", size = 1, linetype = "solid")) + 
  xlab("log2 fold change") + ylab("-log10 q-value") +
  scale_colour_manual(values = c("red", "grey", "steelblue")) + 
  scale_x_continuous(limits = c(-55, 55)) + scale_y_continuous(limits = c(0, 28))

volcDataT <- read.table('deseq2_results_ecic-tissue_LvsR.txt', sep = '\t', header = TRUE)
volcDataT$diffexpr <- "NO"
volcDataT$diffexpr[volcDataT$log2FoldChange > 0 & volcDataT$padj < 0.05] <- "LEAF"
volcDataT$diffexpr[volcDataT$log2FoldChange < 0 & volcDataT$padj < 0.05] <- "ROOT"
volcPlotT <- ggplot(data = volcDataT, aes(x = log2FoldChange, y = -log10(padj), colour = diffexpr)) +	
  geom_point(alpha = 0.4, size = 1.75) + 
  theme_minimal() + labs(title = "Leaf vs Root") +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5), axis.line = element_line(colour = "black", size = 1, linetype = "solid")) + 
  xlab("log2 fold change") + ylab("-log10 q-value") +
  scale_colour_manual(values = c("red", "grey", "steelblue")) + 
  scale_x_continuous(limits = c(-50, 50)) + scale_y_continuous(limits = c(0, 300))

# Age
volcDataYO <- read.table('deseq2_results_ecic-age-YvsO.txt', sep = '\t', header = TRUE)
volcDataYO$diffexpr <- "NO"
volcDataYO$diffexpr[volcDataYO$log2FoldChange > 0 & volcDataYO$padj < 0.05] <- "YOUNG"
volcDataYO$diffexpr[volcDataYO$log2FoldChange < 0 & volcDataYO$padj < 0.05] <- "ADULT"
volcPlotYO <- ggplot(data = volcDataYO, aes(x = log2FoldChange, y = -log10(padj), colour = diffexpr)) +	
  geom_point(alpha = 0.4, size = 1.75) + 
  theme_minimal() + labs(title = "Young vs Adult") +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5), axis.line = element_line(colour = "black", size = 1, linetype = "solid")) + 
  xlab("log2 fold change") + ylab("-log10 q-value") +
  scale_colour_manual(values = c("steelblue", "grey", "red")) + # red = up in Young, blue = up in Old
  scale_x_continuous(limits = c(-25, 25)) + scale_y_continuous(limits = c(0, 80))


volcDataL <- read.table('deseq2_results_ecic-aza_CLYvsALY.txt', sep = '\t', header = TRUE)
volcDataL$diffexpr <- "NO"
volcDataL$diffexpr[volcDataL$log2FoldChange > 0 & volcDataL$padj < 0.05] <- "DOWN"
volcDataL$diffexpr[volcDataL$log2FoldChange < 0 & volcDataL$padj < 0.05] <- "UP"
volcPlotL <- ggplot(data = volcDataL, aes(x = log2FoldChange, y = -log10(padj), colour = diffexpr)) +	
  geom_point(alpha = 0.4, size = 1.75) + 
  theme_minimal() + labs(title = "Control vs Azacytidine (leaf)") +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5), axis.line = element_line(colour = "black", size = 1, linetype = "solid")) + 
  xlab("log2 fold change") + ylab("-log10 q-value") +
  scale_colour_manual(values = c("red", "grey", "steelblue")) + 
  scale_x_continuous(limits = c(-25, 25)) + scale_y_continuous(limits = c(0, 15))

volcDataR <- read.table('deseq2_results_ecic-aza_CRYvsARY.txt', sep = '\t', header = TRUE)
volcDataR$diffexpr <- "NO"
volcDataR$diffexpr[volcDataR$log2FoldChange > 0 & volcDataR$padj < 0.05] <- "DOWN"
volcDataR$diffexpr[volcDataR$log2FoldChange < 0 & volcDataR$padj < 0.05] <- "UP"
volcPlotR <- ggplot(data = volcDataR, aes(x = log2FoldChange, y = -log10(padj), colour = diffexpr)) +	
  geom_point(alpha = 0.4, size = 1.75) + 
  theme_minimal() + labs(title = "Control vs Azacytidine (root)") +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5), axis.line = element_line(colour = "black", size = 1, linetype = "solid")) + 
  xlab("log2 fold change") + ylab("-log10 q-value") +
  scale_colour_manual(values = c("red", "grey", "steelblue")) + 
  scale_x_continuous(limits = c(-30, 30)) + scale_y_continuous(limits = c(0, 30))

volcDataCT <- read.table('deseq2_results_ecic-tissue_CLYvsCRY.txt', sep = '\t', header = TRUE)
volcDataCT$diffexpr <- "NO"
volcDataCT$diffexpr[volcDataCT$log2FoldChange > 0 & volcDataCT$padj < 0.05] <- "LEAF"
volcDataCT$diffexpr[volcDataCT$log2FoldChange < 0 & volcDataCT$padj < 0.05] <- "ROOT"
volcPlotCT <- ggplot(data = volcDataCT, aes(x = log2FoldChange, y = -log10(padj), colour = diffexpr)) +	
  geom_point(alpha = 0.4, size = 1.75) + 
  theme_minimal() + labs(title = "Leaf vs Root (control)") +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5), axis.line = element_line(colour = "black", size = 1, linetype = "solid")) + 
  xlab("log2 fold change") + ylab("-log10 q-value") +
  scale_colour_manual(values = c("red", "grey", "steelblue")) + 
  scale_x_continuous(limits = c(-30, 30)) + scale_y_continuous(limits = c(0, 300))

volcDataAT <- read.table('deseq2_results_ecic-tissue_ALYvsARY.txt', sep = '\t', header = TRUE)
volcDataAT$diffexpr <- "NO"
volcDataAT$diffexpr[volcDataAT$log2FoldChange > 0 & volcDataAT$padj < 0.05] <- "LEAF"
volcDataAT$diffexpr[volcDataAT$log2FoldChange < 0 & volcDataAT$padj < 0.05] <- "ROOT"
volcPlotAT <- ggplot(data = volcDataAT, aes(x = log2FoldChange, y = -log10(padj), colour = diffexpr)) +	
  geom_point(alpha = 0.4, size = 1.75) + 
  theme_minimal() + labs(title = "Leaf vs Root (azacytitdine)") +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5), axis.line = element_line(colour = "black", size = 1, linetype = "solid")) + 
  xlab("log2 fold change") + ylab("-log10 q-value") +
  scale_colour_manual(values = c("red", "grey", "steelblue")) + 
  scale_x_continuous(limits = c(-30, 30)) + scale_y_continuous(limits = c(0, 300))

# Azacytidine within young (same as azacytidine within leaf)
volcDataYCA <- read.table('deseq2_results_ecic-aza_CLYvsALY.txt', sep = '\t', header = TRUE)
volcDataYCA$diffexpr <- "NO"
volcDataYCA$diffexpr[volcDataYCA$log2FoldChange > 0 & volcDataYCA$padj < 0.05] <- "DOWN"
volcDataYCA$diffexpr[volcDataYCA$log2FoldChange < 0 & volcDataYCA$padj < 0.05] <- "UP"
volcPlotYCA <- ggplot(data = volcDataYCA, aes(x = log2FoldChange, y = -log10(padj), colour = diffexpr)) +	
  geom_point(alpha = 0.4, size = 1.75) + 
  theme_minimal() + labs(title = "Control vs Azacytidine (young)") +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5), axis.line = element_line(colour = "black", size = 1, linetype = "solid")) + 
  xlab("log2 fold change") + ylab("-log10 q-value") +
  scale_colour_manual(values = c("red", "grey", "steelblue")) + # red = up in Aza, blue = down in Aza
  scale_x_continuous(limits = c(-25, 25)) + scale_y_continuous(limits = c(0, 4))

# Azacytidine within adult
volcDataOCA <- read.table('deseq2_results_ecic-age-CLOvsALO.txt', sep = '\t', header = TRUE)
volcDataOCA$diffexpr <- "NO"
volcDataOCA$diffexpr[volcDataOCA$log2FoldChange > 0 & volcDataOCA$padj < 0.05] <- "DOWN"
volcDataOCA$diffexpr[volcDataOCA$log2FoldChange < 0 & volcDataOCA$padj < 0.05] <- "UP"
volcPlotOCA <- ggplot(data = volcDataOCA, aes(x = log2FoldChange, y = -log10(padj), colour = diffexpr)) +	
  geom_point(alpha = 0.4, size = 1.75) + 
  theme_minimal() + labs(title = "Control vs Azacytidine (adult)") +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5), axis.line = element_line(colour = "black", size = 1, linetype = "solid")) + 
  xlab("log2 fold change") + ylab("-log10 q-value") +
  scale_colour_manual(values = c("red", "grey", "steelblue")) + # red = up in Aza, blue = down in Aza
  scale_x_continuous(limits = c(-10, 10)) + scale_y_continuous(limits = c(0, 4))

# Age within control
volcDataCYO <- read.table('deseq2_results_ecic-age-CLYvsCLO.txt', sep = '\t', header = TRUE)
volcDataCYO$diffexpr <- "NO"
volcDataCYO$diffexpr[volcDataCYO$log2FoldChange > 0 & volcDataCYO$padj < 0.05] <- "YOUNG"
volcDataCYO$diffexpr[volcDataCYO$log2FoldChange < 0 & volcDataCYO$padj < 0.05] <- "ADULT"
volcPlotCYO <- ggplot(data = volcDataCYO, aes(x = log2FoldChange, y = -log10(padj), colour = diffexpr)) +	
  geom_point(alpha = 0.4, size = 1.75) + 
  theme_minimal() + labs(title = "Young vs Adult (control)") +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5), axis.line = element_line(colour = "black", size = 1, linetype = "solid")) + 
  xlab("log2 fold change") + ylab("-log10 q-value") +
  scale_colour_manual(values = c("steelblue", "grey", "red")) + # red = up in Young, blue = up in Adult
  scale_x_continuous(limits = c(-25, 25)) + scale_y_continuous(limits = c(0, 80))

# Age within azacytidine
volcDataAYO <- read.table('deseq2_results_ecic-age-ALYvsALO.txt', sep = '\t', header = TRUE)
volcDataAYO$diffexpr <- "NO"
volcDataAYO$diffexpr[volcDataAYO$log2FoldChange > 0 & volcDataAYO$padj < 0.05] <- "YOUNG"
volcDataAYO$diffexpr[volcDataAYO$log2FoldChange < 0 & volcDataAYO$padj < 0.05] <- "ADULT"
volcPlotAYO <- ggplot(data = volcDataAYO, aes(x = log2FoldChange, y = -log10(padj), colour = diffexpr)) +	
  geom_point(alpha = 0.4, size = 1.75) + 
  theme_minimal() + labs(title = "Young vs Adult (azacytidine)") +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5), axis.line = element_line(colour = "black", size = 1, linetype = "solid")) + 
  xlab("log2 fold change") + ylab("-log10 q-value") +
  scale_colour_manual(values = c("steelblue", "grey", "red")) + # red = up in Young, blue = up in Adult
  scale_x_continuous(limits = c(-25, 25)) + scale_y_continuous(limits = c(0, 80))

library(gridExtra)
# all general effects
grid.arrange(volcPlotA, volcPlotT, volcPlotYO, ncol = 3, nrow = 1)
# tissue differences
grid.arrange(volcPlotL, volcPlotR, volcPlotCT, volcPlotAT, ncol = 2, nrow = 2)
# age differences
grid.arrange(volcPlotYCA, volcPlotOCA, volcPlotCYO, volcPlotAYO, ncol = 2, nrow = 2)
# within group differences
grid.arrange(volcPlotR, volcPlotL, volcPlotOCA, ncol = 3, nrow = 1)


#### FISHER'S EXACT TEST ####
#phyper(q,m,n,k)
#q = size of overlap - 1
#m = number of DEGs in experiment comparison 1
#n = total number of genes - m (70973 - m)
#k = number of DEGs in comparison experiment 2

## General effects ##
# Intersection CvsA and LvsR: 100
length(intersect(as.vector(rownames(data.CvsA[which(data.CvsA$padj <= 0.05),])), as.vector(rownames(data.LvsR[which(data.LvsR$padj <= 0.05),]))))
# Intersection CvsA and YvsO: 85
length(intersect(as.vector(rownames(data.CvsA[which(data.CvsA$padj <= 0.05),])), as.vector(rownames(data.YvsO[which(data.YvsO$padj <= 0.05),]))))
# Intersection LvsR and YvsO: 29379
length(intersect(as.vector(rownames(data.LvsR[which(data.LvsR$padj <= 0.05),])), as.vector(rownames(data.YvsO[which(data.YvsO$padj <= 0.05),]))))

# CvsA (164) and LvsR (54424)
phyper(99, 164, 67930, 54424, lower.tail = FALSE, log.p = FALSE) # P-value = 1
# CvsA (164) and YvsO (37102)
phyper(84, 164, 67930, 37102, lower.tail = FALSE, log.p = FALSE) # P-value = 0.777
# LvsR (54424) and YvsO (37102)
phyper(29378, 54424, 13670, 37102, lower.tail = FALSE, log.p = FALSE) # P-value = 0.999

