#!/usr/bin/env nextflow

// Specify DSL2
nextflow.preview.dsl = 2

// Process definition
process deseq2 {
    publishDir "${params.outdir}/deseq2",
        mode: "copy", overwrite: true

    input:
        path(count_data)

    output:
        path(results)

    script:

    deseq2_command = "Rscript $baseDir/deseq2.R -c $count_data -m ${params.metadata}"
    
    // SHELL
    """
    ${deseq2_command}
    """
}