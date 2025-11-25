#### E. CICUTARIUM AZACITIDINE: COUNTS PREPARATION ####

setwd("/path/to/project/directory/results/WGCNA")

library('WGCNA')
source('wgcna2igraph.R')
library('dynamicTreeCut')
library(stats)
library(stringr)

# FORMATTING NORMALIZED COUNTS

normFilteredCounts <- read.csv("../deseq2/ecic-all_normcounts_deseq2.txt", sep = "\t", header = T)

data.Expr <- t(log2(normFilteredCounts + 1))
data.Gene <- rownames(normFilteredCounts)
data.Sample <- colnames(normFilteredCounts) 
colnames(data.Expr) <- data.Gene
rownames(data.Expr) <- data.Sample

# POWER OF THE SOFT THRESHOLD VALUE

beta1 <- 6 #for 20-30 samples, 8 is recommended for unsigned networks	
k.dat.Expr <- softConnectivity(data.Expr, power = beta1) - 1	
scaleFreePlot(k.dat.Expr, main = paste("data Expression, power = ", beta1), truncated = F)
dev.copy(pdf, file = "connectivity_distribution.pdf")
dev.off()

#### GENERATING DENDROGRAM ####

kCut <- (dim(data.Expr)[[2]])
kRank <- rank(-k.dat.Expr)
vardat.Expr <- apply(data.Expr, 2, var) 
restk <- kRank <= kCut & vardat.Expr > 0 

ADJdat.Expr <- adjacency(datExpr = data.Expr[,restk], power = beta1, corOptions = "use = 'p', method = 'spearman'")	#Use corOptions = "use = 'p', method = 'spearman'" in further analysis since is less a less sensitive method to outliers compared to Pearson (default method)
dissTOMdat.Expr <- TOMdist(ADJdat.Expr)
hierTOMdat.Expr <- hclust(as.dist(dissTOMdat.Expr), method = "average")	#Use as input for cutreeDynamic()

# STATIC METHOD MODULES (DISREGARD)

colorhdat.Expr <- cutreeStaticColor(hierTOMdat.Expr, cutHeight = 0.9, minSize = 50)	#defaults 0.9 and minSize = 50

#### MODULE DETECTION: HYBRID METHOD ####

# DETECT MODULES: check how parameters detectCutHeight and deepSplit shape the modules

# data.net1 <- blockwiseModules(datExpr = data.Expr[,restk], maxBlockSize = 70973, networkType = "unsigned",
#                              power = beta1, detectCutHeight = 0.98, deepSplit = 0, minModuleSize = 50, 
#                              saveTOMs = FALSE, verbose = F) #maxBlockSize = 70973
# data.net2 <- blockwiseModules(datExpr = data.Expr[,restk], maxBlockSize = 70973, networkType = "unsigned",
#                              power = beta1, detectCutHeight = 0.98, deepSplit = 4, minModuleSize = 50, 
#                              saveTOMs = FALSE, verbose = F) #maxBlockSize = 70973
# data.net3 <- blockwiseModules(datExpr = data.Expr[,restk], maxBlockSize = 70973, networkType = "unsigned",
#                              power = beta1, detectCutHeight = 0.9, deepSplit = 0, minModuleSize = 50, 
#                              saveTOMs = FALSE, verbose = F) #maxBlockSize = 70973
data.net4 <- blockwiseModules(datExpr = data.Expr[,restk], maxBlockSize = 70973, networkType = "unsigned",
                             power = beta1, detectCutHeight = 0.9, deepSplit = 4, minModuleSize = 50, 
                             saveTOMs = FALSE, verbose = F) #maxBlockSize = 70973

# table(data.net$colors)
# genesByModule1 <- data.frame(cbind(data.Gene[restk], data.net1$colors))
# genesByModule2 <- data.frame(cbind(data.Gene[restk], data.net2$colors))
# genesByModule3 <- data.frame(cbind(data.Gene[restk], data.net3$colors))
genesByModule4 <- data.frame(cbind(data.Gene[restk], data.net4$colors))
#write.table(genesByModule, file = "WGCNA_genesByModule_defaultOptions.txt", sep = "\t")
# write.table(genesByModule1, file = "WGCNA_genesByModule_d098s0m050.txt", sep = "\t")
# write.table(genesByModule2, file = "WGCNA_genesByModule_d098s4m050.txt", sep = "\t")
# write.table(genesByModule3, file = "WGCNA_genesByModule_d090s0m050.txt", sep = "\t")
write.table(genesByModule4, file = "WGCNA_genesByModule_d090s4m050.txt", sep = "\t")

# PLOT MODULES IN THE DENDROGRAM

# hybTOMdat.Expr <- data.net$colors
# hybTOMdat1.Expr <- data.net1$colors
# hybTOMdat2.Expr <- data.net2$colors
# hybTOMdat3.Expr <- data.net3$colors
hybTOMdat4.Expr <- data.net4$colors

# par(mfrow=c(6,1),mar=c(2,4,1,1))
# plot(hierTOMdat.Expr, main = "Network Dendrogram", labels = F, xlab = "", sub = "")
# plotColorUnderTree(hierTOMdat.Expr, colors = data.frame(module = colorhdat.Expr)) #static
# plotColorUnderTree(hierTOMdat.Expr, colors = data.frame(module = hybTOMdat1.Expr)) #hybrid
# plotColorUnderTree(hierTOMdat.Expr, colors = data.frame(module = hybTOMdat2.Expr)) #hybrid
# plotColorUnderTree(hierTOMdat.Expr, colors = data.frame(module = hybTOMdat3.Expr)) #hybrid
# plotColorUnderTree(hierTOMdat.Expr, colors = data.frame(module = hybTOMdat4.Expr)) #hybrid
# dev.copy(pdf, "WGCNA_moduleMethodComparison_all.pdf")
# dev.off()

data.net <- data.net4
hybTOMdat.Expr <- hybTOMdat4.Expr # Module set (N) to be used
genesByModule <- genesByModule4

par(mfrow=c(3,1),mar=c(2,4,1,1))
plot(hierTOMdat.Expr, main = "Network Dendrogram", labels = F, xlab = "", sub = "")
plotColorUnderTree(hierTOMdat.Expr, colors = data.frame(module = colorhdat.Expr)) #static
plotColorUnderTree(hierTOMdat.Expr, colors = data.frame(module = hybTOMdat.Expr)) #hybrid
dev.copy(pdf, "WGCNA_moduleMethodComparison.pdf")
dev.off()

#### DETECT HUB GENES ####

# Top hub genes per module (use the hub_score() function after getting the igraph object)
# topHubs <- as.matrix(chooseTopHubInEachModule(data.Expr[,restk], data.net$colors, power = beta1, type = "unsigned"))
# write.table(topHubs, file = "WGCNA_topHubs.txt", sep = "\t")

blockColors <- rownames(sort(table(data.net$colors), decreasing = T))

graph <- wgcna2igraph(net = data.net, datExpr = data.Expr[, restk],
                      modules2plot = blockColors,
                      colors2plot = blockColors,
                      kME.threshold = 0.5, adjacency.threshold = 0.5, adj.power = beta1,
                      verbose = T, node.size = 1, frame.color = "black", node.color = "red",
                      edge.alpha = .5, edge.width = 1)

# plot(graph)

# Compute Kleinberg's hub centrality scores
hub.genes <- hub_score(graph, scale = T)
hub.genes$vector
write.table(hub.genes$vector, file = "WGCNA_hub_scores.txt", sep = "\t")


#### EIGENVALUE STATISTICS ####

data.net$MEs$azacytidine <- factor(str_sub(as.character(rownames(data.Expr)), -5, -5))
data.net$MEs$tissue <- factor(str_sub(as.character(rownames(data.Expr)), -3, -3))
data.net$MEs$age <- factor(str_sub(as.character(rownames(data.Expr)), -1, -1))
write.table(data.net$MEs, file = "MEs_values_for_stats.txt", sep = "\t", quote = F)

# Use when data.net is not loaded in R
# data.net <- list()
# data.net$MEs <- read.table("MEs_values_for_stats.txt", sep = '\t', header = T, row.names = 1)

# PCA
rownames(data.net$MEs) <- rownames(data.Expr)
MEpca <- prcomp(data.net$MEs[1:12])
s <- summary(MEpca) # Proportion of variance: PC1 = 52.8%, PC2 = 19.4%
# Color vector: white = CLY, red = AYL, orange = CRY, brown = ARY, green = CLO, darkgreen = ALO
col.group<- c("white", "orange", "red", "brown", 
             "white", "orange", "red", "brown", 
             "white", "red", "brown", "white", 
             "orange", "red", "brown", "white", 
             "orange", "red", "brown", "white", 
             "orange", "red", "brown", "green", 
             "darkgreen", "green", "darkgreen", "green", 
             "darkgreen", "darkgreen", "green", "darkgreen", 
             "green", "darkgreen", "green")

plot(MEpca$x[,1], MEpca$x[,2], 
     xlab = paste("PCA 1 (", round(s$importance[2]*100, 1), "%)", sep = ""), 
     ylab = paste("PCA 2 (", round(s$importance[5]*100, 1), "%)", sep = ""), 
     pch = 21, col = "black", bg = col.group, cex = 1, las = 1, asp = 1)
abline(v = 0, lty = 2, col = "grey50")
abline(h = 0, lty = 2, col = "grey50")
l.x <- MEpca$rotation[,1] * 1.25
l.y <- MEpca$rotation[,2] * 1.25
arrows(x0 = 0, x1 = l.x, y0 = 0, y1 = l.y, col = "red", length = 0.15, lwd = 1.5)
l.pos <- l.y
lo <- which(l.y < 0)
hi <- which(l.y > 0)
l.pos <- replace(l.pos, lo, "1")
l.pos <- replace(l.pos, hi, "3")
text(l.x, l.y, labels = row.names(MEpca$rotation), col = "red", pos = l.pos)

MEcolors <- colnames(data.net$MEs[1:16]) # range from 1 to total number of modules
MEcounter <- 1

# KRUSKAL-WALLIS & DUNN TESTS

while (MEcounter <= length(MEcolors)){
  Azacytidine.kw <- kruskal.test(as.formula(paste0("data.net$MEs$", MEcolors[MEcounter], " ~ data.net$MEs$azacytidine")))
  Tissue.kw <- kruskal.test(as.formula(paste0("data.net$MEs$", MEcolors[MEcounter], " ~ data.net$MEs$tissue")))
  Age.kw <- kruskal.test(as.formula(paste0("data.net$MEs$", MEcolors[MEcounter], " ~ data.net$MEs$age")))
  MEvalues <- as.data.frame(data.net$MEs[MEcounter])
  my_title <- MEcolors[MEcounter]
  my_list <- list(my_title, Azacytidine.kw, Tissue.kw, Age.kw)
  capture.output(my_list, file = "ME_KWtest_results.txt", append = T)
  MEcounter <- MEcounter + 1
}

# BOXPLOT FOR ALL MODULES

library(ggplot2)
library(ggpubr)
library(RColorBrewer)


MEcolors <- colnames(data.net$MEs)[1:13]
data.net$MEs$age.tissue <- paste(data.net$MEs$age, data.net$MEs$tissue, sep = ".")
MEcounter <- 1

pdf(file = "MEs_boxplots_19042024.pdf")
while (MEcounter <= length(data.net$MEs)){
  bxplt <- ggboxplot(data.net$MEs, x = "age.tissue", y = MEcolors[MEcounter], color = "azacytidine", 
                     repel = T, font.label = list(size = 14, face = "plain"), 
                     add = "jitter", shape = "azacytidine",
                     ggtheme = theme_gray())
  print(bxplt)
  MEcounter <- MEcounter + 1
}
dev.off()

# CORRELATION BETWEEN MODULES

MEcor <- cor(as.data.frame(data.net$MEs[1:12])) # range from 1 to total number of modules

library(pheatmap)
library(RColorBrewer)

pheatmap(MEcor, color = colorRampPalette(brewer.pal(n = 7, name = "RdBu"))(100))

#### GENE SIGNIFICANCE (-log(FDR)) & GENE MODULE MEMBERSHIP (cor(gene, ME)) ####

# FUNCTION TO GENERATE A MEAN + SE SUMMARY
data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm = TRUE),
      se = sd(x[[col]], na.rm = TRUE)/sqrt(length(x[[col]])))
  }
  data_sum <- ddply(data, groupnames, .fun = summary_func,
                    varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}

# DATA INPUT
exprGenes <- read.table("ecic-all_normcounts_deseq2.txt", sep = "\t", header = T, row.names = 1)
listModules <- read.table("WGCNA_genesByModule_d090s4m050.txt", sep = "\t", header = T, row.names = 1)
listModules$X1 <- NULL
colnames(listModules) <- "module"
exprGenes$module <- listModules$module 

# CALCULATE GENE SIGNIFICANCE PER MODULE
listPvalues <- read.csv("WGCNA_data.csv", sep = ";", header = T, row.names = 1) # compile FDRs from DESeq2
exprGenes <- exprGenes[order(row.names(exprGenes)), ]
listPvalues <- listPvalues[order(row.names(listPvalues)), ]
exprGenes$geneSig.aza <- -log(as.numeric(listPvalues$FDR.aza) + 0.001, 2) # FDR from control vs azacytidine contrast
exprGenes$geneSig.tis <- -log(as.numeric(listPvalues$FDR.tissue) + 0.001, 2) # FDR from leaf vs root contrast
exprGenes$geneSig.age <- -log(as.numeric(listPvalues$FDR.age) + 0.001, 2) # FDR from young vs adult contrast
exprGenes$geneSig.roo <- -log(as.numeric(listPvalues$FDR.root) + 0.001, 2) # FDR from control vs azacytidine within roots
exprGenes$geneSig.lea <- -log(as.numeric(listPvalues$FDR.leaf) + 0.001, 2) # FDR from control vs azacytidine within leaves
exprGenes$geneSig.adu <- -log(as.numeric(listPvalues$FDR.adult) + 0.001, 2) # FDR from control vs azacytidine within adults

blockColors <- names(sort(table(exprGenes$module), decreasing = T))
blockColors[1:2] <- c("grey", "turquoise") # Use to swap M1 and singlets when N genes in singlets is greater than in M1

# azacytidine
gsdata.aza <- as.data.frame(cbind(exprGenes$module, as.numeric(exprGenes$geneSig.aza)))
rownames(gsdata.aza) <- rownames(exprGenes)
colnames(gsdata.aza) <- c("module", "geneSig.aza")
gsdata.aza$module <- factor(gsdata.aza$module, 
                            levels = blockColors)
levels(gsdata.aza$module) <- c("singlets", "M1", "M2", "M3", "M4", "M5", "M6", 
                               "M7", "M8", "M9", "M10", "M11", "M12") # Adjust module vector
gsdata.aza$geneSig.aza <- as.numeric(gsdata.aza$geneSig.aza)
# tissue
gsdata.tis <- as.data.frame(cbind(exprGenes$module, as.numeric(exprGenes$geneSig.tis)))
rownames(gsdata.tis) <- rownames(exprGenes)
colnames(gsdata.tis) <- c("module", "geneSig.tis")
gsdata.tis$module <- factor(gsdata.tis$module, 
                            levels = blockColors)
levels(gsdata.tis$module) <- c("singlets", "M1", "M2", "M3", "M4", "M5", "M6", 
                               "M7", "M8", "M9", "M10", "M11", "M12") # Adjust module vector
gsdata.tis$geneSig.tis <- as.numeric(gsdata.tis$geneSig.tis)
# age
gsdata.age <- as.data.frame(cbind(exprGenes$module, as.numeric(exprGenes$geneSig.age)))
rownames(gsdata.age) <- rownames(exprGenes)
colnames(gsdata.age) <- c("module", "geneSig.age")
gsdata.age$module <- factor(gsdata.age$module, 
                            levels = blockColors)
levels(gsdata.age$module) <- c("singlets", "M1", "M2", "M3", "M4", "M5", "M6", 
                               "M7", "M8", "M9", "M10", "M11", "M12") # Adjust module vector
gsdata.age$geneSig.age <- as.numeric(gsdata.age$geneSig.age)
# roo
gsdata.roo <- as.data.frame(cbind(exprGenes$module, as.numeric(exprGenes$geneSig.roo)))
rownames(gsdata.roo) <- rownames(exprGenes)
colnames(gsdata.roo) <- c("module", "geneSig.roo")
gsdata.roo$module <- factor(gsdata.roo$module, 
                            levels = blockColors)
levels(gsdata.roo$module) <- c("singlets", "M1", "M2", "M3", "M4", "M5", "M6", 
                               "M7", "M8", "M9", "M10", "M11", "M12") # Adjust module vector
gsdata.roo$geneSig.roo <- as.numeric(gsdata.roo$geneSig.roo)
# lea
gsdata.lea <- as.data.frame(cbind(exprGenes$module, as.numeric(exprGenes$geneSig.lea)))
rownames(gsdata.lea) <- rownames(exprGenes)
colnames(gsdata.lea) <- c("module", "geneSig.lea")
gsdata.lea$module <- factor(gsdata.lea$module, 
                            levels = blockColors)
levels(gsdata.lea$module) <- c("singlets", "M1", "M2", "M3", "M4", "M5", "M6", 
                               "M7", "M8", "M9", "M10", "M11", "M12") # Adjust module vector
gsdata.lea$geneSig.lea <- as.numeric(gsdata.lea$geneSig.lea)
# adu
gsdata.adu <- as.data.frame(cbind(exprGenes$module, as.numeric(exprGenes$geneSig.adu)))
rownames(gsdata.adu) <- rownames(exprGenes)
colnames(gsdata.adu) <- c("module", "geneSig.adu")
gsdata.adu$module <- factor(gsdata.adu$module, 
                            levels = blockColors)
levels(gsdata.adu$module) <- c("singlets", "M1", "M2", "M3", "M4", "M5", "M6", 
                               "M7", "M8", "M9", "M10", "M11", "M12") # Adjust module vector
gsdata.adu$geneSig.adu <- as.numeric(gsdata.adu$geneSig.adu)

# drop singlets by choosing rows 2 to 13 in each summary 
gsdata.aza.dsum <- data_summary(gsdata.aza, varname = "geneSig.aza", groupnames = "module")[2:13,]
gsdata.tis.dsum <- data_summary(gsdata.tis, varname = "geneSig.tis", groupnames = "module")[2:13,]
gsdata.age.dsum <- data_summary(gsdata.age, varname = "geneSig.age", groupnames = "module")[2:13,]

gsdata.roo.dsum <- data_summary(gsdata.roo, varname = "geneSig.roo", groupnames = "module")[2:13,]
gsdata.lea.dsum <- data_summary(gsdata.lea, varname = "geneSig.lea", groupnames = "module")[2:13,]
gsdata.adu.dsum <- data_summary(gsdata.adu, varname = "geneSig.adu", groupnames = "module")[2:13,]

library(ggplot2)
library(RColorBrewer)
# work in a custom palete where module color is representative
# azacytidine
gsdata.aza.plot <- ggplot(gsdata.aza.dsum, aes(x = module, y = geneSig.aza, fill = module)) + 
                   geom_bar(stat = "identity", color = "black", 
                            position = position_dodge()) +
                   geom_errorbar(aes(ymin = geneSig.aza - se, ymax = geneSig.aza + se), width=.2,
                                 position = position_dodge(.9)) +
                   scale_fill_manual("Module", values = blockColors[2:13])
# tissue
gsdata.tis.plot <- ggplot(gsdata.tis.dsum, aes(x = module, y = geneSig.tis, fill = module)) + 
  geom_bar(stat = "identity", color = "black", 
           position = position_dodge()) +
  geom_errorbar(aes(ymin = geneSig.tis - se, ymax = geneSig.tis + se), width=.2,
                position = position_dodge(.9)) +
  scale_fill_manual("Module", values = blockColors[2:13])
# age
gsdata.age.plot <- ggplot(gsdata.age.dsum, aes(x = module, y = geneSig.age, fill = module)) + 
  geom_bar(stat = "identity", color = "black", 
           position = position_dodge()) +
  geom_errorbar(aes(ymin = geneSig.age - se, ymax = geneSig.age + se), width=.2,
                position = position_dodge(.9)) +
  scale_fill_manual("Module", values = blockColors[2:13])

library(gridExtra)
grid.arrange(gsdata.aza.plot, gsdata.tis.plot, gsdata.age.plot, ncol = 1, nrow = 3)

# root
gsdata.roo.plot <- ggplot(gsdata.roo.dsum, aes(x = module, y = geneSig.roo, fill = module)) + 
  geom_bar(stat = "identity", color = "black", 
           position = position_dodge()) +
  geom_errorbar(aes(ymin = geneSig.roo - se, ymax = geneSig.roo + se), width=.2,
                position = position_dodge(.9)) +
  scale_fill_manual("Module", values = blockColors[2:13])
# leaf
gsdata.lea.plot <- ggplot(gsdata.lea.dsum, aes(x = module, y = geneSig.lea, fill = module)) + 
  geom_bar(stat = "identity", color = "black", 
           position = position_dodge()) +
  geom_errorbar(aes(ymin = geneSig.lea - se, ymax = geneSig.lea + se), width=.2,
                position = position_dodge(.9)) +
  scale_fill_manual("Module", values = blockColors[2:13])
# adult
gsdata.adu.plot <- ggplot(gsdata.adu.dsum, aes(x = module, y = geneSig.adu, fill = module)) + 
  geom_bar(stat = "identity", color = "black", 
           position = position_dodge()) +
  geom_errorbar(aes(ymin = geneSig.adu - se, ymax = geneSig.adu + se), width=.2,
                position = position_dodge(.9)) +
  scale_fill_manual("Module", values = blockColors[2:13])

library(gridExtra)
grid.arrange(gsdata.roo.plot, gsdata.lea.plot, gsdata.adu.plot, ncol = 1, nrow = 3)

# CALCULATE GENE MODULE MEMBERSHIP

moduleEigen <- read.table("MEs_values_for_stats.txt", sep = '\t', header = T) # exported from data.net$MEs, load only if not already

## Plot with GS vs GMM plots for azacytidine
par(mfcol = c(2,2))

#  Module 6 (red)
gsdata.m6 <- gsdata.aza[which(gsdata.aza$module == "M6"), ]
exprGenes.m6 <- exprGenes[which(exprGenes$module == "red"), ]

i <- 1
gModMem <- vector()
while (i <= length(rownames(exprGenes.m6))) {
  gModMem[i] <- cor(as.numeric(exprGenes.m6[i, c(1:35)]), moduleEigen$MEred)
  i <- i + 1
}

plot(as.numeric(gsdata.m6$geneSig.aza), as.numeric(abs(gModMem)),
     xlab = "M6 gene significance", ylab = "M6 |gene module membership|", col = "red")

#  Module 10 (purple)
# gsdata.m10 <- gsdata.aza[which(gsdata.aza$module == "M10"), ]
# exprGenes.m10 <- exprGenes[which(exprGenes$module == "purple"), ]
# 
# i <- 1
# gModMem <- vector()
# while (i <= length(rownames(exprGenes.m10))) {
#   gModMem[i] <- cor(as.numeric(exprGenes.m10[i, c(1:35)]), moduleEigen$MEpurple)
#   i <- i + 1
# }
# 
# plot(as.numeric(gsdata.m10$geneSig.aza), as.numeric(abs(gModMem)),
#      xlab = "M10 gene significance", ylab = "M10 |gene module membership|", col = "purple")

#  Module 9 (magenta)
gsdata.m9 <- gsdata.aza[which(gsdata.aza$module == "M9"), ]
exprGenes.m9 <- exprGenes[which(exprGenes$module == "magenta"), ]

i <- 1
gModMem <- vector()
while (i <= length(rownames(exprGenes.m9))) {
  gModMem[i] <- cor(as.numeric(exprGenes.m9[i, c(1:35)]), moduleEigen$MEmagenta)
  i <- i + 1
}

plot(as.numeric(gsdata.m9$geneSig.aza), as.numeric(abs(gModMem)),
     xlab = "M9 gene significance", ylab = "M9 |gene module membership|", col = "magenta")

#  Module 11 (greenyellow)
gsdata.m11 <- gsdata.aza[which(gsdata.aza$module == "M11"), ]
exprGenes.m11 <- exprGenes[which(exprGenes$module == "greenyellow"), ]

i <- 1
gModMem <- vector()
while (i <= length(rownames(exprGenes.m11))) {
  gModMem[i] <- cor(as.numeric(exprGenes.m11[i, c(1:35)]), moduleEigen$MEgreenyellow)
  i <- i + 1
}

plot(as.numeric(gsdata.m11$geneSig.aza), as.numeric(abs(gModMem)),
     xlab = "M11 gene significance", ylab = "M11 |gene module membership|", col = "greenyellow")

#  Module 12 (tan)
gsdata.m12 <- gsdata.aza[which(gsdata.aza$module == "M12"), ]
exprGenes.m12 <- exprGenes[which(exprGenes$module == "tan"), ]

i <- 1
gModMem <- vector()
while (i <= length(rownames(exprGenes.m12))) {
  gModMem[i] <- cor(as.numeric(exprGenes.m12[i, c(1:35)]), moduleEigen$MEtan)
  i <- i + 1
}

plot(as.numeric(gsdata.m12$geneSig.aza), as.numeric(abs(gModMem)),
     xlab = "M12 gene significance", ylab = "M12 |gene module membership|", col = "tan")

dev.off()

## Other GS vs GMM plots
#  Module 8 (pink) with within-root FDR
gsdata.m8 <- gsdata.roo[which(gsdata.roo$module == "M8"), ]
exprGenes.m8 <- exprGenes[which(exprGenes$module == "pink"), ]

i <- 1
gModMem <- vector()
while (i <= length(rownames(exprGenes.m8))) {
  gModMem[i] <- cor(as.numeric(exprGenes.m8[i, c(1:35)]), moduleEigen$MEpink)
  i <- i + 1
}

plot(as.numeric(gsdata.m8$geneSig.roo), as.numeric(abs(gModMem)),
     xlab = "M8 gene significance", ylab = "M8 |gene module membership|", col = "pink")

#  Module 12 (tan) with within.leaf FDR
gsdata.m12 <- gsdata.lea[which(gsdata.lea$module == "M12"), ]
exprGenes.m12 <- exprGenes[which(exprGenes$module == "tan"), ]

i <- 1
gModMem <- vector()
while (i <= length(rownames(exprGenes.m12))) {
  gModMem[i] <- cor(as.numeric(exprGenes.m12[i, c(1:35)]), moduleEigen$MEtan)
  i <- i + 1
}

plot(as.numeric(gsdata.m12$geneSig.lea), as.numeric(abs(gModMem)),
     xlab = "M12 gene significance", ylab = "M12 |gene module membership|", col = "tan")

#  Module 9 (magenta) with within.adult FDR
gsdata.m9 <- gsdata.lea[which(gsdata.lea$module == "M9"), ]
exprGenes.m9 <- exprGenes[which(exprGenes$module == "magenta"), ]

i <- 1
gModMem <- vector()
while (i <= length(rownames(exprGenes.m9))) {
  gModMem[i] <- cor(as.numeric(exprGenes.m9[i, c(1:35)]), moduleEigen$MEmagenta)
  i <- i + 1
}

plot(as.numeric(gsdata.m9$geneSig.lea), as.numeric(abs(gModMem)),
     xlab = "M9 gene significance", ylab = "M9 |gene module membership|", col = "magenta")
