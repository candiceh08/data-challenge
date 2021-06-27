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

include {deseq2} from 'rna-seq/deseq2/main.nf' 

workflow rnaseq-da {
    main:
        // Normalise the count data

        // 
}