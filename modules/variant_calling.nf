process VARIANT_CALLING {
    container 'broadinstitute/gatk:latest'
    publishDir "${params.outdir}/variants", mode: 'copy'

    input:
    tuple val(sample), path(bam)
    path "*"
    path reference
    path "*"
    path "*"

    output:
    tuple val(sample), path("variants.vcf.gz")

    script:
    """
    gatk HaplotypeCaller \
        -R ${reference} \
        -I ${bam} \
        -O variants.vcf.gz
    """
}
