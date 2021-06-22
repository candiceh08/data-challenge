#!/usr/bin/env nextflow

// Specify DSL2
nextflow.enable.dsl=2

// Log 
log.info("Starting downstream analysis of rna-seq and atac-seq data")

/*------------------------------------------------------------------------------------*/
/* Define params
--------------------------------------------------------------------------------------*/



/*------------------------------------------------------------------------------------*/
/* Module inclusions
--------------------------------------------------------------------------------------*/

include {normalisation} from 'rna-seq/normalisation/main.nf' 
include {dea} from 'rna-seq/diff-expression-analysis/main.nf'
include{uca} from 'rna-seq/unsupervised-clustering-analyses/main.nf'



workflow {
  
}