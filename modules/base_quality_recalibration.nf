process BASE_QUALITY_RECALIBRATION {
    tag "${sample}"
    label 'low_mem'
    container 'broadinstitute/gatk:latest'
    publishDir "${params.outdir}/recalibrated", mode: 'copy'

    input:
    tuple val(sample), path(bam), path(bam_idx)
    tuple path(reference), path("*"), path("*"), path("*")
    path known_sites

    output:
    tuple val(sample), path("recalibrated.bam")

    script:
    """
    if [ -z "${known_sites}" ]; then
        echo "WARNING: No known-sites found. Skipping BaseRecalibrator."
        cp ${bam} recalibrated.bam
    else
        gatk IndexFeatureFile \
            -I ${known_sites} \
            --java-options "-Xmx${task.memory.toGiga()}g -XX:+UseParallelGC -XX:ParallelGCThreads=${task.cpus}"

        gatk BaseRecalibrator \
            -I ${bam} \
            -R ${reference} \
            --known-sites ${known_sites} \
            -O recalibration.table \
            --java-options "-Xmx${task.memory.toGiga()}g -XX:+UseParallelGC -XX:ParallelGCThreads=${task.cpus}"

        gatk ApplyBQSR \
            -I ${bam} \
            -R ${reference} \
            --bqsr-recal-file recalibration.table \
            -O recalibrated.bam \
            --java-options "-Xmx${task.memory.toGiga()}g -XX:+UseParallelGC -XX:ParallelGCThreads=${task.cpus}"
    fi
    """
}
