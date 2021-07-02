#!/usr/bin/env nextflow

// Specify DSL2
nextflow.preview.dsl = 2

// Process definition
process deseq2 {
    publishDir "${params.outdir}/deseq2",
        mode: "copy", 
    saveAs: {filename ->
        if (filename.indexOf(".txt") > 0) "deseq2_PlotsAndFiles/files/$filename"
        else if (filename == "PCAplot.png") "deseq2_PlotsAndFiles/$filename"
        else if (filename.indexOf("_MAplot.png") > 0) "deseq2_PlotsAndFiles/MAplots/$filename"
        else "$filename"
    }

    container "deseq2"

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