# Data Challenge

Genomics and Regulatory Systems Unit (Luscombe Unit), OIST

### Introduction

Here you will find the data from an RNA-Seq and ATAC-Seq experiments. Both data sets have been generated using the same experimental design. There is a treatment and control group each containing three replicates making a total of six samples per experiment. The data is accessible via Dropbox links below. The data files are defined as follows (all files are tab delimited text files):

#### RNA-Seq Data

rnaseq_design.txt		Sample ids and corresponding condition labels.
rnaseq_gene_counts.txt	Raw (not normalised) gene-level read counts for each sample.
rnaseq_annotation.txt		Gene level annotation.

#### ATAC-Seq Data

atacseq_design.txt		Sample ids and corresponding condition labels.
atacseq_peak_counts.txt	Raw (not normalised) ATAC-Seq peak level counts for each sample.
atacseq_peaks.bed		A bed file defining the peak loci

All sequence data were aligned to the human genome reference hg38. 
Data can be downloaded from this here https://www.dropbox.com/s/xij5qe3xl1uul4m/data.zip

### The Challenge

The treatment is thought to activate a transcriptional program via remodelling of the chromatin architecture. Without using a premade pipeline to perform these analyses, please do the following:

    1) Identify genes that may be regulated in this fashion.
    2) Identify the possible transcriptional programs involved.
    3) Present candidate transcription factors that may be responsible for the underlying regulation.
    4) Implement your analysis such that it could be reproduced by someone else.

Please send a link to any source code before the interview. For the interview, please prepare a 30-minute presentation detailing your exploration of the data, your analysis approach, your findings, and the implementation of your solution. 

