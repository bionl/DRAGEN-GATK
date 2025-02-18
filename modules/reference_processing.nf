process REFERENCE_DICT {
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
