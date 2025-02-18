include { BUILD_HASH_TABLE } from './modules/build_hash_table'
include { ALIGN_READS } from './modules/align_reads'
include { CONVERT_SAM_TO_BAM_AND_SORT } from './modules/convert_sam_to_bam'
include { REFERENCE_DICT ; REFERENCE_INDEXING } from './modules/reference_processing'
include { BASE_QUALITY_RECALIBRATION } from './modules/base_quality_recalibration'
include { BAM_INDEXING } from './modules/bam_indexing'
include { VARIANT_CALLING } from './modules/variant_calling'

workflow {
    /**
     * Load input sample sheet and reference genome
     * - `params.input`: Path to a CSV file containing sample information
     * - `params.reference`: Path to the reference genome
     * - `params.known_sites`: (Optional) Path to known variant sites for base quality recalibration
     */
    samples = Channel.fromPath(params.input).splitCsv(header: true)
    reference = file(params.reference)
    known_sites = params.known_sites ? file(params.known_sites) : null

    /**
     * Reference genome preprocessing
     * - `BUILD_HASH_TABLE`: Creates a hash table for reference genome indexing
     * - `REFERENCE_DICT`: Generates a sequence dictionary for the reference
     * - `REFERENCE_INDEXING`: Creates an index file for the reference genome
     */
    hash_table = BUILD_HASH_TABLE(reference)
    ref_dict = REFERENCE_DICT(reference)
    ref_index = REFERENCE_INDEXING(reference)

    /**
     * Read alignment step
     * - `ALIGN_READS`: Aligns sequencing reads to the reference genome
     * - Requires reference genome and precomputed hash table
     */
    aligned_reads = ALIGN_READS(samples, reference, hash_table)

    /**
     * Conversion and sorting
     * - `CONVERT_SAM_TO_BAM_AND_SORT`: Converts SAM to BAM format and sorts the reads
     */
    sorted_bam = CONVERT_SAM_TO_BAM_AND_SORT(aligned_reads)

    /**
     * Base Quality Score Recalibration (BQSR)
     * - Only performed if known sites are provided
     * - `BASE_QUALITY_RECALIBRATION`: Adjusts base quality scores to correct sequencing errors
     */
    if (known_sites) {
        recalibrated_bam = BASE_QUALITY_RECALIBRATION(sorted_bam, reference, hash_table, ref_dict, ref_index, known_sites)
    }
    else {
        recalibrated_bam = sorted_bam
    }

    /**
     * BAM indexing step
     * - `BAM_INDEXING`: Generates an index file for the final BAM file
     */
    indexed_bam = BAM_INDEXING(recalibrated_bam)

    /**
     * Variant calling
     * - `VARIANT_CALLING`: Identifies variants from the aligned and processed reads
     * - Requires indexed BAM file and reference genome
     */
    VARIANT_CALLING(recalibrated_bam, indexed_bam, reference, ref_dict, ref_index)
}
