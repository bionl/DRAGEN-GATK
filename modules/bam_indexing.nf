process BAM_INDEXING {
    container 'staphb/samtools:latest'

    input:
    tuple val(sample), path(bam)

    output:
    path "*.bai"

    script:
    """
    samtools index ${bam}
    """
}
