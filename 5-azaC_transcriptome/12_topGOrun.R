#### TopGO: GO TERM ENRICHMENT ANALYSIS ####

library(topGO)

setwd("/path/to/project/directory/results/topGO")

# Load a background GO term set
geneID2GO <- readMappings(file = "GOterm_universe_ecic-RNAseq.txt")

length(geneID2GO) # number of genes in the GO universe

# Define colMap helper function
colMap <- function(x) {
  .col <- rep(rev(heat.colors(length(unique(x)))), time = table(x))
  return(.col[match(1:length(x), order(x))])
}

# Generate a list of files to analyze
degList <- list.files(recursive = F, pattern = "ecic-RNAseq*")

# degList <- c("m10root.txt", "m9.txt", "m12.txt") # modules from WGCNA

#### topGO LOOP ####
# Run in $PATH with gene lists with adjusted p-value
for (i in 1:length(degList)) {
  GOlist <- t(read.csv(degList[[i]], sep = "\t", header = F))
  GOlist.v <- as.numeric(c(GOlist[2,]))
  names(GOlist.v) <- GOlist[1,]
  # Generate factor TRUE/FALSE for the genes in gene list present in whole gene set
  geneNames <- names(geneID2GO)
  myInterestingGenes <- GOlist[1,]
  GOlist.list <- factor(as.integer(geneNames %in% myInterestingGenes))
  names(GOlist.list) <- geneNames
  ## BIOLOGICAL PROCESS ##
  # Create topGOdata object
  GOlist_BP <- new("topGOdata", 
                   description = "", 
                   ontology = "BP", # which GO term family to analyze
                   allGenes = GOlist.list, # input gene list
                   nodeSize = 5, # filter nodes (GO terms) with less than 5 seq counts
                   annot = annFUN.gene2GO, # function for extracting reference annotation from GO mapping file
                   gene2GO = geneID2GO) # mapping file
  # Classic Fisher exact test
  GOlist_BP.fet <- runTest(GOlist_BP, algorithm = "classic", statistic = "fisher")
  # Weighted Fisher exact test
  GOlist_BP.weight <- runTest(GOlist_BP, algorithm = "weight", statistic = "fisher")
  # Compile and save the results
  GOlist_BP.allRes <- GenTable(GOlist_BP, classicFisher = GOlist_BP.fet, weightFisher = GOlist_BP.weight,
                               orderBy = "weightFisher", ranksOf = "weightFisher", topNodes = 11524) # BP: 11524 nodes
  GOlist_BP.allRes$classicFisher[GOlist_BP.allRes$classicFisher == "< 1e-30"] <- "1e-30"
  GOlist_BP.classicFDR <- p.adjust(GOlist_BP.allRes$classicFisher, method = "fdr")
  GOlist_BP.allRes$classicFDR <- GOlist_BP.classicFDR
  GOlist_BP.allRes$weightFisher[GOlist_BP.allRes$weightFisher == "< 1e-30"] <- "1e-30"
  GOlist_BP.weightFDR <- p.adjust(GOlist_BP.allRes$weightFisher, method = "fdr")
  GOlist_BP.allRes$weightFDR <- GOlist_BP.weightFDR
  write.table(GOlist_BP.allRes, 
              file = paste0(degList[[i]], "_topGO_results_BP.txt"), 
              sep = '\t', 
              quote = F)
  ## MOLECULAR FUNCTION ##
  # Create topGOdata object
  GOlist_MF <- new("topGOdata", 
                   description = "", 
                   ontology = "MF", # which GO term family to analyze
                   allGenes = GOlist.list, # input gene list
                   nodeSize = 5, # filter nodes (GO terms) with less than 5 seq counts
                   annot = annFUN.gene2GO, # function for extracting reference annotation from GO mapping file
                   gene2GO = geneID2GO) # mapping file
  # Classic Fisher exact test
  GOlist_MF.fet <- runTest(GOlist_MF, algorithm = "classic", statistic = "fisher")
  # Weighted Fisher exact test
  GOlist_MF.weight <- runTest(GOlist_MF, algorithm = "weight", statistic = "fisher")
  # Compile and save the results
  GOlist_MF.allRes <- GenTable(GOlist_MF, classicFisher = GOlist_MF.fet, weightFisher = GOlist_MF.weight,
                               orderBy = "weightFisher", ranksOf = "weightFisher", topNodes = 3317) # MF: 3317 nodes
  GOlist_MF.allRes$classicFisher[GOlist_MF.allRes$classicFisher == "< 1e-30"] <- "1e-30"
  GOlist_MF.classicFDR <- p.adjust(GOlist_MF.allRes$classicFisher, method = "fdr")
  GOlist_MF.allRes$classicFDR <- GOlist_MF.classicFDR
  GOlist_MF.allRes$weightFisher[GOlist_MF.allRes$weightFisher == "< 1e-30"] <- "1e-30"
  GOlist_MF.weightFDR <- p.adjust(GOlist_MF.allRes$weightFisher, method = "fdr")
  GOlist_MF.allRes$weightFDR <- GOlist_MF.weightFDR
  write.table(GOlist_MF.allRes, 
              file = paste0(degList[[i]], "_topGO_results_MF.txt"), 
              sep = '\t', 
              quote = F)
  ## CELLULAR COMPONENT ##
  # Create topGOdata object
  GOlist_CC <- new("topGOdata", 
                   description = "", 
                   ontology = "CC", # which GO term family to analyze
                   allGenes = GOlist.list, # input gene list
                   nodeSize = 5, # filter nodes (GO terms) with less than 5 seq counts
                   annot = annFUN.gene2GO, # function for extracting reference annotation from GO mapping file
                   gene2GO = geneID2GO) # mapping file
  # Classic Fisher exact test
  GOlist_CC.fet <- runTest(GOlist_CC, algorithm = "classic", statistic = "fisher")
  # Weighted Fisher exact test
  GOlist_CC.weight <- runTest(GOlist_CC, algorithm = "weight", statistic = "fisher")
  # Compile and save the results
  GOlist_CC.allRes <- GenTable(GOlist_CC, classicFisher = GOlist_CC.fet, weightFisher = GOlist_CC.weight,
                               orderBy = "weightFisher", ranksOf = "weightFisher", topNodes = 1617) # CC: 1617 nodes
  GOlist_CC.allRes$classicFisher[GOlist_CC.allRes$classicFisher == "< 1e-30"] <- "1e-30"
  GOlist_CC.classicFDR <- p.adjust(GOlist_CC.allRes$classicFisher, method = "fdr")
  GOlist_CC.allRes$classicFDR <- GOlist_CC.classicFDR
  GOlist_CC.allRes$weightFisher[GOlist_CC.allRes$weightFisher == "< 1e-30"] <- "1e-30"
  GOlist_CC.weightFDR <- p.adjust(GOlist_CC.allRes$weightFisher, method = "fdr")
  GOlist_CC.allRes$weightFDR <- GOlist_CC.weightFDR
  write.table(GOlist_CC.allRes, 
              file = paste0(degList[[i]], "_topGO_results_CC.txt"), 
              sep = '\t', 
              quote = F)
}

#### REVIGO: PARAMETER DESCRIPTION ####

# Run ReViGO analysis at: http://revigo.irb.hr/

# Analyze lists with a large number of enriched GO terms (> 10), from biological process

# use Small (0.5) word similarity threshold, and UniProtKB as reference database

# Integrate ReViGO and topGO data, and select top ten enriched GO terms from each list 
# (choose GO terms classified as non redundant by ReViGO) 

#### GO TERM ENRICHMENT PLOTS ####

# Prepare a compilation of enriched GO terms from several contrasts in Excel and import to R
goterms <- read.csv("GOtermsSummary.csv", sep = ';', header = T)
goterms <- goterms[which(goterms$Ontology == "CC"), ] # Ontology options: "BP", "MF", "CC"

## enriched GO term lists ##

goterms.root <- goterms[which(goterms$comparison == "CRvsAR_AR"), ]
goterms.adult <- goterms[which(goterms$comparison == "COvsAO_AO"), ]
goterms.module6 <- goterms[which(goterms$comparison == "module6"), ]

require(ggplot2)
library(scales)
library(forcats)

goplot.root <- ggplot(goterms.root,
                      aes(x = fct_reorder(Term, weightFDR, .desc = TRUE), y = -log10(weightFDR), size = log2(Significant/Expected), fill = -log10(weightFDR))) +
  expand_limits(y = 1) +
  geom_point(shape = 21) +
  scale_size(range = c(2.5,12.5)) +
  scale_fill_continuous(low = 'yellow', high = 'red') +
  
  xlab('') + ylab('Enrichment score') +
  labs(title = 'Enriched GO terms (root-azacytidine)') +
  
  theme_bw(base_size = 24) +
  theme(
    legend.position = 'right',
    legend.background = element_rect(),
    plot.title = element_text(angle = 0, size = 16, face = 'bold', vjust = 1),
    plot.subtitle = element_text(angle = 0, size = 14, face = 'bold', vjust = 1),
    plot.caption = element_text(angle = 0, size = 12, face = 'bold', vjust = 1),
    
    axis.text.x = element_text(angle = 0, size = 12, face = 'bold', hjust = 1.10),
    axis.text.y = element_text(angle = 0, size = 12, face = 'bold', vjust = 0.5),
    axis.title = element_text(size = 12, face = 'bold'),
    axis.title.x = element_text(size = 12, face = 'bold'),
    axis.title.y = element_text(size = 12, face = 'bold'),
    axis.line = element_line(colour = 'black'),
    
    #Legend
    legend.key = element_blank(), # removes the border
    legend.key.size = unit(1, "cm"), # Sets overall area/size of the legend
    legend.text = element_text(size = 14, face = "bold"), # Text size
    title = element_text(size = 14, face = "bold")) +
  
  coord_flip()

goplot.adult <- ggplot(goterms.adult,
                       aes(x = fct_reorder(Term, weightFDR, .desc = TRUE), y = -log10(weightFDR), size = log2(Significant/Expected), fill = -log10(weightFDR))) +
  expand_limits(y = 1) +
  geom_point(shape = 21) +
  scale_size(range = c(2.5,12.5)) +
  scale_fill_continuous(low = 'yellow', high = 'red') +
  
  xlab('') + ylab('Enrichment score') +
  labs(title = 'Enriched GO terms (adult-azacytidine)') +
  
  theme_bw(base_size = 24) +
  theme(
    legend.position = 'right',
    legend.background = element_rect(),
    plot.title = element_text(angle = 0, size = 16, face = 'bold', vjust = 1),
    plot.subtitle = element_text(angle = 0, size = 14, face = 'bold', vjust = 1),
    plot.caption = element_text(angle = 0, size = 12, face = 'bold', vjust = 1),
    
    axis.text.x = element_text(angle = 0, size = 12, face = 'bold', hjust = 1.10),
    axis.text.y = element_text(angle = 0, size = 12, face = 'bold', vjust = 0.5),
    axis.title = element_text(size = 12, face = 'bold'),
    axis.title.x = element_text(size = 12, face = 'bold'),
    axis.title.y = element_text(size = 12, face = 'bold'),
    axis.line = element_line(colour = 'black'),
    
    #Legend
    legend.key = element_blank(), # removes the border
    legend.key.size = unit(1, "cm"), # Sets overall area/size of the legend
    legend.text = element_text(size = 14, face = "bold"), # Text size
    title = element_text(size = 14, face = "bold")) +
  
  coord_flip()

goplot.module6 <- ggplot(goterms.module6,
                         aes(x = fct_reorder(Term, weightFDR, .desc = TRUE), y = -log10(weightFDR), size = log2(Significant/Expected), fill = -log10(weightFDR))) +
  expand_limits(y = 1) +
  geom_point(shape = 21) +
  scale_size(range = c(2.5,12.5)) +
  scale_fill_continuous(low = 'yellow', high = 'red') +
  
  xlab('') + ylab('Enrichment score') +
  labs(title = 'Enriched GO terms (module6-azacytidine)') +
  
  theme_bw(base_size = 24) +
  theme(
    legend.position = 'right',
    legend.background = element_rect(),
    plot.title = element_text(angle = 0, size = 16, face = 'bold', vjust = 1),
    plot.subtitle = element_text(angle = 0, size = 14, face = 'bold', vjust = 1),
    plot.caption = element_text(angle = 0, size = 12, face = 'bold', vjust = 1),
    
    axis.text.x = element_text(angle = 0, size = 12, face = 'bold', hjust = 1.10),
    axis.text.y = element_text(angle = 0, size = 12, face = 'bold', vjust = 0.5),
    axis.title = element_text(size = 12, face = 'bold'),
    axis.title.x = element_text(size = 12, face = 'bold'),
    axis.title.y = element_text(size = 12, face = 'bold'),
    axis.line = element_line(colour = 'black'),
    
    #Legend
    legend.key = element_blank(), # removes the border
    legend.key.size = unit(1, "cm"), # Sets overall area/size of the legend
    legend.text = element_text(size = 14, face = "bold"), # Text size
    title = element_text(size = 14, face = "bold")) +
  
  coord_flip()


library(gridExtra)
goplot.root2 <- ggplotGrob(goplot.root)
goplot.adult2 <- ggplotGrob(goplot.adult)
goplot.module62 <- ggplotGrob(goplot.module6)
grid.arrange(goplot.root2, goplot.adult2, goplot.module62, nrow = 3, ncol = 1) # 12 in x 37 in

