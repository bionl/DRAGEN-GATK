process BUILD_HASH_TABLE {
    label 'big_mem'
    container 'gambalab/dragmap:latest'

    input:
    path reference

    output:
    path "*"

    script:
    """
    dragen-os  \
        --build-hash-table true \
        --ht-reference ${reference} \
        --output-directory ./
    """
}
