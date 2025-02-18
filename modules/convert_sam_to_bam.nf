process CONVERT_SAM_TO_BAM_AND_SORT {
    container 'staphb/samtools:latest'

    input:
    tuple val(sample), path(x)

    output:
    tuple val(sample), path("sorted.bam")

    script:
    """
    # Ensure output format is BAM
    if [ -f "${sample}_aligned.sam" ]; then
        samtools view -bS ${sample}_aligned.sam > ${sample}_aligned.bam
    fi

    # Sort & Index BAM
    samtools sort -o sorted.bam ${sample}_aligned.bam
    samtools index sorted.bam
    """
}
