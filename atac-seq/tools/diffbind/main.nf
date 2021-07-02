#!/usr/bin/env nextflow

// Specify DSL2
nextflow.preview.dsl = 2

// Process definition
process diffbind {
    publishDir "${params.outdir}/deseq2", mode: "copy", overwrite = true

    container "biowardrobe2/diffbind"

    input:
        path(count_data)

    output:
        path'dds_results.txt', emit: dds_results
        path'dds_summary.txt', emit: dds_summary
        path'dds_significant_genes.txt', emit: dds_significant_genes
        path'PCA_plot.pdf', emit: pca_plot

    script:

    deseq2_command = "Rscript $baseDir/deseq2.R -c $count_data -m ${params.metadata} -o ${params.outdir} "
    
       // SHELL
    """
    ${deseq2_command}
    """
}