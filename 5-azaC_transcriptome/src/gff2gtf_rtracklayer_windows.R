# R script to convert NCBI's gff to gtf using rtracklayer package.
# Call from command line on biocluster2 using:
# $ srun --pty /bin/bash # to get an interactive session
# $ module load  R/3.5.0-IGB-gcc-4.9.4
# $ cd /path_to_folder_with_NCBI_gff (uncompressed or gzipped)
# $ Rscript gff2gtf_rtracklayer.R existing.gff outputFileName.gtf putEGin_gene_id

# if specify putEGin_gene_id, the script will pull out the Entrez GeneIDs from
# the Dbxref column and put them into the attribute column under the standard "gene_id"
# Also, if no information in Dbxref column, use info in product column

#### DO NOT RUN WITH RSCRIPT COMMAND!!!
#### RUN COPYING AND PASTING THE CODE ON A R PROMPT

setwd("C:\\Users\\hp\\Documents\\Ruben\\Postdoc\\2018_CameronLabPostdoc\\Genomics")
args <- c("GCF_000188095.3_BIMP_2.2_genomic.gff","GCF_000188095.3_BIMP_2.2_genomic.gtf","putEGin_gene_id") # B. impatiens
#args <- c("GCF_000214255.1_Bter_1.0_genomic.gff","GCF_000214255.1_Bter_1.0_genomic.gtf","putEGin_gene_id") # B. terrestris
gff0 <- rtracklayer::import(args[1])
if (length(args) == 3) {
  if(args[3] == "putEGin_gene_id") {
    temp <- lapply(gff0$Dbxref, function(x) x[grep("GeneID", x)])
    #If didn't have GeneID, switch to NA
    temp2 <- sapply(temp, length)
      if(any(temp2 == 0))
        temp[temp2 == 0] <- NA
    #If have more than one GeneID, keep first one
      if(any(temp2 > 1))
        temp[temp2 > 1] <- lapply(temp[temp2 > 1], function(x) x[1])
    #Put GeneID in attribute named gene_id
      gff0$gene_id <- unlist(temp)
    #If no GeneID, put in value from product attribute
      if(sum(is.na(gff0$gene_id)) > 0)
        gff0$gene_id[is.na(gff0$gene_id)] <- gff0$product[is.na(gff0$gene_id)]
   }
}
rtracklayer::export(gff0, args[2], format = "gtf")

