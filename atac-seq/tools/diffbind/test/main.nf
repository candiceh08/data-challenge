#!/usr/bin/env nextflow

// Define DSL2
nextflow.preview.dsl=2

// Log
log.info ("Starting tests for diffbind...")

/*------------------------------------------------------------------------------------*/
/* Define params
--------------------------------------------------------------------------------------*/

params.metadata = "$baseDir/data/atacseq_design.txt"

/*------------------------------------------------------------------------------------*/
/* Module inclusions
--------------------------------------------------------------------------------------*/

include {diffbind} from '../main.nf'

/*------------------------------------------------------------------------------------*/
/* Define input channels
--------------------------------------------------------------------------------------*/

testcounts = "$baseDir/data/atacseq_peak_counts.txt"

//Define test data input channels
Channel
    .from(testcounts)
    .set{ch_test_counts}

/*------------------------------------------------------------------------------------*/
/* Run tests
--------------------------------------------------------------------------------------*/

workflow {
    // Run deseq2
    deseq2(ch_test_counts)
}