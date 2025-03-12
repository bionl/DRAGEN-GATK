process CONVERT_SAM_TO_BAM_AND_SORT {
    tag "${sample}"
    label "med_mem"
    container 'staphb/samtools:latest'

    input:
    tuple val(sample), path(x)

    output:
    tuple val(sample), path("${sample}.sorted.bam"), path("${sample}.sorted.bam.bai")

    script:
    """
    # Ensure output format is BAM
    if [ -f "${sample}_aligned.sam" ]; then
        samtools view -bS ${sample}_aligned.sam > ${sample}_aligned.bam
    fi

    # Sort & Index BAM
    samtools sort -m ${task.memory.toMega()}M -@ ${task.cpus} -o ${sample}.sorted.bam ${sample}_aligned.bam
    samtools index ${sample}.sorted.bam
    """
}
