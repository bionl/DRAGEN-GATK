process BAM_INDEXING {
    tag "${sample}"
    container 'staphb/samtools:latest'

    input:
    tuple val(sample), path(bam)

    output:
    tuple val(sample), path("*.bai")

    script:
    """
    samtools index ${bam}
    """
}
