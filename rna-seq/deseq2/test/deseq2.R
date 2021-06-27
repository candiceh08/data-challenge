#!/usr/bin/env Rscript

#-------------------------------------------------------------------------------------------------------------------------
#- R script that runs differential gene expression using the DESeq2 bioconductor package
#-------------------------------------------------------------------------------------------------------------------------

# TO-DO: 
# create two options: one for metadata path, one for count data path - done
# add all the R packages to the docker image - done 

# supress warnings
#options(warn=-1)

#suppressPackageStartupMessages(library(optparse))
#suppressPackageStartupMessages(library(ggplot2))
#suppressPackageStartupMessages(library(pheatmap))
#suppressPackageStartupMessages(library(mclust))
#suppressPackageStartupMessages(library(factoextra))
#suppressPackageStartupMessages(library(cowplot))
#suppressPackageStartupMessages(library(gridExtra))
#suppressPackageStartupMessages(library(dplyr))
#suppressPackageStartupMessages(library(stringr))
#suppressPackageStartupMessages(library(openxlsx))
#suppressPackageStartupMessages(library(RColorBrewer))
#suppressPackageStartupMessages(library(colorspace))
#suppressPackageStartupMessages(library(data.table))

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

opt = parse_args(parser)

counts = opt$count_data_file
metadata = opt$metadata_file

#-------------------------------------------------------------------------------------------------------------------------
#- Get metadata
#-------------------------------------------------------------------------------------------------------------------------
colData <- read.table(metadata, sep="\t",header=T)
samples<-colData[,1]
colnames<-colnames(colData)
colData<-as.data.frame(colData[,-c(1)])
rownames(colData)<-samples
colnames(colData)<-colnames[-1]

#------------------------------------------------------------------------------------------------------------------------
#- Create a DESeqDataSet object with DESeqDataSetFromMatrix function
#------------------------------------------------------------------------------------------------------------------------

countData <- read.table(counts)
dds <- DESeqDataSetFromMatrix(countData = countData,
                                colData = colData,
                                design = eval(parse(text=paste0("~ ", design))))

#------------------------------------------------------------------------------------------------------------------------
#- Filtering steps
#------------------------------------------------------------------------------------------------------------------------

#- Remove genes with <= 5 counts in all samples 

dds <- dds[ rowSums(counts(dds)) > 5, ]

#------------------------------------------------------------------------------------------------------------------------
#- Run DGE
#------------------------------------------------------------------------------------------------------------------------

dds <- DESeq(dds)