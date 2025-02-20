process BUILD_HASH_TABLE {
    label 'big_mem'
    container 'gambalab/dragmap:latest'

    input:
    tuple val(filename), path(reference)

    output:
    tuple val(filename), path("*")

    script:
    """
    dragen-os  \
        --build-hash-table true \
        --ht-reference ${reference} \
        --output-directory ./ \
        --num-threads ${task.cpus}
    """
}
