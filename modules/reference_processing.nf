process REFERENCE_DICT {
    label "low_mem"
    container 'broadinstitute/gatk:latest'

    input:
    tuple val(filename), path(reference)

    output:
    tuple val(filename), path("*.dict")

    script:
    """
    gatk CreateSequenceDictionary \
        -R ${reference}
    """
}

process REFERENCE_INDEXING {
    label "low_mem"
    container 'staphb/samtools:latest'

    input:
    tuple val(filename), path(reference)

    output:
    tuple val(filename), path("*.fai")

    script:
    """
    samtools faidx ${reference}
    """
}
