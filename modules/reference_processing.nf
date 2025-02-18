process REFERENCE_DICT {
    label "low_mem"
    container 'broadinstitute/gatk:latest'

    input:
    path reference

    output:
    path "*.dict"

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
    path reference

    output:
    path "*.fai"

    script:
    """
    samtools faidx ${reference}
    """
}
