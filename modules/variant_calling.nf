process VARIANT_CALLING {
    tag "${sample}"
    label 'low_mem'
    container 'broadinstitute/gatk:latest'
    publishDir "${params.outdir}/variants", mode: 'copy'

    input:
    tuple val(sample), path(bam), path(bam_idx)
    tuple path(reference), path("*"), path("*"), path("*")

    output:
    tuple val(sample), path("variants.vcf.gz")

    script:
    """
    gatk \
        --java-options "-Xmx${task.memory.toGiga()}g -XX:+UseParallelGC -XX:ParallelGCThreads=${task.cpus}" \
        HaplotypeCaller \
        -R ${reference} \
        -I ${bam} \
        -O variants.vcf.gz \
        --native-pair-hmm-threads ${task.cpus}
    """
}
