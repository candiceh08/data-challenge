#!/usr/bin/env Rscript

#-------------------------------------------------------------------------------------------------------------------------
#- R script that runs differential gene expression using the DESeq2 bioconductor package
#-------------------------------------------------------------------------------------------------------------------------

#----------------------------------
#- Load / install required packages
#----------------------------------

library(optparse)
library(stringr)
library(BiocParallel)
library(chromVAR)
library(magrittr)
library(GenomicRanges)
library(DESeq2)
library(ggplot2)
library(ggbeeswarm)
library(dplyr)
library(RColorBrewer)
library(pheatmap)
library(AnnotationDbi)
library(org.Hs.eg.db)

#library(lattice)

# Creating options

parser <- OptionParser()

parser <- add_option(parser, 
                     opt_str = c("-c", "--count_data_file"), 
                     type = "character",
                     dest = 'count_data_file',
                     help="Path to the count data file which contains the non-normalized gene-level counts (required)."
                    )

parser <- add_option(parser, 
                     opt_str = c("-m", "--metadata"), 
                     type = "character",
                     dest = 'metadata_file',
                     help="Path to the metadata data file (required)."
                    )
                
parser <- add_option(parser,
                    opt_str = c("-o", "--outdir"), 
                     type = "character",
                     dest = 'output_dir',
                     help="Output directory where results are saved."
                    )

opt = parse_args(parser)

counts = opt$count_data_file
metadata = opt$metadata_file
outdir = opt$output_dir

#-------------------------------------------------------------------------------------------------------------------------
#- Get metadata
#-------------------------------------------------------------------------------------------------------------------------
colData <- read.table(metadata, sep="\t", header=T)
#samples <- colData[,1] #makes a vector of all the values in the first column of the table
#columnnames <-colnames(colData)
#colData<-as.data.frame(colData[,-c(1)])
#rownames(colData) <- samples
#colnames(colData) <- columnnames[-1]

#------------------------------------------------------------------------------------------------------------------------
#- Create a DESeqDataSet object with DESeqDataSetFromMatrix function
#------------------------------------------------------------------------------------------------------------------------

# Create a table for counts
countData <- read.table(counts, header = T,sep = "\t",check.names = FALSE)
#geneID <- countData$featureid
#countData <- select(countData, -featureid)
#rownames(countData) <- geneID

dds <- DESeqDataSetFromMatrix(countData = countData,
                                colData = colData,
                                design = ~condition, tidy = TRUE)

#------------------------------------------------------------------------------------------------------------------------
#- Filtering steps
#------------------------------------------------------------------------------------------------------------------------

#- Remove genes with <= 5 counts in all samples 

dds <- dds[ rowSums(counts(dds)) > 5, ]

#------------------------------------------------------------------------------------------------------------------------
#- Run DGE and get results table for most significant genes
#------------------------------------------------------------------------------------------------------------------------

dds <- DESeq(dds)
res <- results(dds, tidy=TRUE) # Creates a data frame containing the results of DGE
write.table(res,"dds_results.txt", append=FALSE, sep="\t", row.names = TRUE, col.names = TRUE) # Transforms the data frame into a readable text file, which is passed as one of the results of the process deseq2

summary <- summary(res, 0.01, tidy=TRUE) # Gives a summary of the contrasts effects between control and treatment group. Sets the adjusted p-value cutoff at 0.1
write.table(summary,"dds_summary.txt", sep = "\t",row.names = FALSE)

resSig <- res[ which(res$padj < 0.01 ), ] #Filters the reads that have an adjusted p-value lower than 0.01 (1%). padj corresponds to the fraction of the genes that would be false positives if one considered all the genes with a higher p-value than this specific gene to be significant
resSig <- resSig[ order( resSig$log2FoldChange ), ]
write.table(resSig,"dds_significant_genes.txt", append=FALSE, sep="\t", row.names = TRUE, col.names = TRUE)

#------------------------------------------------------------------------------------------------------------------------
#- Add column for the Entrez IDs -> KEGG pathway analysis
#------------------------------------------------------------------------------------------------------------------------

resSig$symbol = mapIds(org.Hs.eg.db,
                     keys=row.names(res), 
                     column="SYMBOL",
                     keytype="ENSEMBL",
                     multiVals="first")
resSig$entrez = mapIds(org.Hs.eg.db,
                     keys=row.names(res), 
                     column="ENTREZID",
                     keytype="ENSEMBL",
                     multiVals="first")
resSig$name =   mapIds(org.Hs.eg.db,
                     keys=row.names(res), 
                     column="GENENAME",
                     keytype="ENSEMBL",
                     multiVals="first")

#------------------------------------------------------------------------------------------------------------------------
#- Plot counts for top 6 differentially expressed genes
#------------------------------------------------------------------------------------------------------------------------
#par(mfrow=c(2,3))

#plotCounts(dds, gene="ENSG00000283186", intgroup="condition")
#plotCounts(dds, gene="ENSG00000283266", intgroup="condition")
#plotCounts(dds, gene="ENSG00000157654", intgroup="condition")
#plotCounts(dds, gene="ENSG00000166349", intgroup="condition")
#plotCounts(dds, gene="ENSG00000267472", intgroup="condition")
#plotCounts(dds, gene="ENSG00000163884", intgroup="condition")

#------------------------------------------------------------------------------------------------------------------------
#- PCA
#------------------------------------------------------------------------------------------------------------------------

#rlogdata <- rlog(dds, blind=FALSE)
#pcaData <- plotPCA(rlogdata,intgroup=c("condition"),ntop = dim(rlogdata)[1], returnData=TRUE)
#percentVar <- round(100*attr(pcaData, "percentVar"))
#pca <- ggplot(pcaData, aes(PC1, PC2, color=colData[,1])) +
#  geom_point(size=3) +
 # xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  #ylab(paste0("PC2: ",percentVar[2], "% variance")) +
  #theme(legend.title = element_blank()) +
  #coord_fixed()
#ggsave(plot = pca, filename = "PCA_plot.pdf", device = "pdf", dpi = 300)
#ggsave(plot = pca, filename = "PCA_plot.svg", device = "svg", dpi = 150)