process ALIGN_READS {
    tag "${sample}"
    label 'mid_mem'
    container 'gambalab/dragmap:latest'
    publishDir "${params.outdir}/aligned_reads", mode: 'copy'

    input:
    tuple val(sample), path(fastq_1), path(fastq_2)
    tuple path(reference), path("*"), path("*"), path("*")

    output:
    tuple val(sample), path("${sample}_aligned*")

    script:
    """
    dragen-os \
        -r ./ \
        -1 ${fastq_1} \
        -2 ${fastq_2} \
        --output-directory ./ \
        --output-file-prefix ${sample}_aligned \
        --num-threads ${task.cpus}
    """
}
