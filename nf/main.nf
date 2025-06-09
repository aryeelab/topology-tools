#!/usr/bin/env nextflow

/*
========================================================================================
    Micro-C Pipeline - Nextflow DSL2 Implementation
========================================================================================
    A comprehensive Nextflow pipeline for processing Micro-C sequencing data.
    
    Based on the Python implementation in microc_pipeline.py with the following steps:
    1. BWA alignment with pairtools parsing
    2. Merge, deduplication, and splitting of pairs
    3. Cooler contact matrix generation (raw and balanced mcool files)
    4. Juicer HiC file generation (with error handling)
    5. QC metrics extraction
    
    Author: Converted from Python pipeline
    Version: 1.0.0
========================================================================================
*/

nextflow.enable.dsl = 2

/*
========================================================================================
    PARAMETER VALIDATION
========================================================================================
*/

// Check required parameters
if (!params.sample_id) {
    error "Parameter 'sample_id' is required"
}
if (!params.fastq_r1) {
    error "Parameter 'fastq_r1' is required"
}
if (!params.fastq_r2) {
    error "Parameter 'fastq_r2' is required"
}
if (!params.reference_bwa_idx) {
    error "Parameter 'reference_bwa_idx' is required"
}
if (!params.chrom_sizes) {
    error "Parameter 'chrom_sizes' is required"
}

// Set default parameters
params.output_dir = "./results"
params.resolution = 1000
params.bwa_cores = 2
params.mapq = 20
params.help = false

/*
========================================================================================
    HELP MESSAGE
========================================================================================
*/

def helpMessage() {
    log.info"""
    ========================================================================================
                            Micro-C Pipeline - Nextflow DSL2
    ========================================================================================
    
    Usage:
        nextflow run nf/main.nf --sample_id SAMPLE --fastq_r1 R1.fq.gz --fastq_r2 R2.fq.gz \\
                                --reference_bwa_idx genome.tgz --chrom_sizes genome.sizes
    
    Required Parameters:
        --sample_id             Sample identifier
        --fastq_r1              R1 FASTQ file(s) (comma-separated for multiple files)
        --fastq_r2              R2 FASTQ file(s) (comma-separated for multiple files)
        --reference_bwa_idx     BWA index tar.gz file
        --chrom_sizes           Chromosome sizes file
    
    Optional Parameters:
        --output_dir            Output directory (default: ./results)
        --resolution            Contact matrix resolution in bp (default: 1000)
        --bwa_cores             Number of BWA cores (default: 2)
        --mapq                  Minimum mapping quality (default: 20)
    
    Profiles:
        -profile local          Run locally (default)
        -profile cluster        Run on cluster with SLURM
        -profile cloud          Run on cloud with appropriate settings
        -profile conda          Use conda environments
        -profile docker         Use Docker containers
        -profile singularity    Use Singularity containers
    
    Examples:
        # Basic usage
        nextflow run nf/main.nf --sample_id sample1 \\
                                --fastq_r1 reads_R1.fq.gz \\
                                --fastq_r2 reads_R2.fq.gz \\
                                --reference_bwa_idx genome.tgz \\
                                --chrom_sizes genome.sizes

        # Multiple FASTQ pairs
        nextflow run nf/main.nf --sample_id sample1 \\
                                --fastq_r1 "reads1_R1.fq.gz,reads2_R1.fq.gz" \\
                                --fastq_r2 "reads1_R2.fq.gz,reads2_R2.fq.gz" \\
                                --reference_bwa_idx genome.tgz \\
                                --chrom_sizes genome.sizes

        # High-performance run
        nextflow run nf/main.nf --sample_id large_sample \\
                                --fastq_r1 reads_R1.fq.gz \\
                                --fastq_r2 reads_R2.fq.gz \\
                                --reference_bwa_idx genome.tgz \\
                                --chrom_sizes genome.sizes \\
                                --bwa_cores 16 --resolution 5000 \\
                                -profile cluster
    
    ========================================================================================
    """.stripIndent()
}

// Show help message if requested
if (params.help) {
    helpMessage()
    exit 0
}

/*
========================================================================================
    IMPORT MODULES
========================================================================================
*/

include { EXTRACT_BWA_INDEX } from './modules/extract_bwa_index'
include { BWA_ALIGN_PAIRTOOLS } from './modules/bwa_align_pairtools'
include { MERGE_DEDUP_SPLIT } from './modules/merge_dedup_split'
include { COOLER_PROCESS } from './modules/cooler_process'
include { JUICER_HIC } from './modules/juicer_hic'
include { QC_METRICS } from './modules/qc_metrics'

/*
========================================================================================
    MAIN WORKFLOW
========================================================================================
*/

workflow {
    
    // Print pipeline information
    log.info """
    ========================================================================================
                            Micro-C Pipeline - Nextflow DSL2
    ========================================================================================
    Sample ID       : ${params.sample_id}
    FASTQ R1        : ${params.fastq_r1}
    FASTQ R2        : ${params.fastq_r2}
    BWA Index       : ${params.reference_bwa_idx}
    Chrom Sizes     : ${params.chrom_sizes}
    Output Dir      : ${params.output_dir}
    Resolution      : ${params.resolution}
    BWA Cores       : ${params.bwa_cores}
    Min MAPQ        : ${params.mapq}
    ========================================================================================
    """.stripIndent()
    
    // Create input channels and pair R1/R2 files correctly
    fastq_r1_list = params.fastq_r1.split(',').collect()
    fastq_r2_list = params.fastq_r2.split(',').collect()

    // Ensure we have the same number of R1 and R2 files
    if (fastq_r1_list.size() != fastq_r2_list.size()) {
        error "Number of R1 files (${fastq_r1_list.size()}) does not match number of R2 files (${fastq_r2_list.size()})"
    }

    // Create paired channel by zipping R1 and R2 files
    fastq_pairs_ch = Channel.fromList(
        fastq_r1_list.indices.collect { i ->
            def r1 = file(fastq_r1_list[i])
            def r2 = file(fastq_r2_list[i])
            def chunk_id = r1.baseName.replaceAll(/\.fq\.gz$|\.fastq\.gz$/, '')
            [chunk_id, r1, r2]
        }
    )

    // Create other input channels
    bwa_index_ch = Channel.fromPath(params.reference_bwa_idx)
    chrom_sizes_ch = Channel.fromPath(params.chrom_sizes)

    // Extract BWA index for each chunk (needed for parallel processing)
    EXTRACT_BWA_INDEX(
        fastq_pairs_ch.map { chunk_id, r1, r2 -> chunk_id },
        bwa_index_ch.first()
    )
    
    // Align reads and process with pairtools
    BWA_ALIGN_PAIRTOOLS(
        fastq_pairs_ch,
        EXTRACT_BWA_INDEX.out.genome_index,
        chrom_sizes_ch.first()
    )
    
    // Collect all pairsam files for merging
    pairsam_files_ch = BWA_ALIGN_PAIRTOOLS.out.pairsam
        .map { chunk_id, pairsam_file -> pairsam_file }
        .collect()
        .map { pairsams -> [params.sample_id, pairsams] }
    
    // Merge, deduplicate, and split pairs
    MERGE_DEDUP_SPLIT(
        pairsam_files_ch
    )
    
    // Process with cooler
    COOLER_PROCESS(
        MERGE_DEDUP_SPLIT.out.mapped_pairs,
        chrom_sizes_ch.first()
    )

    // Generate Juicer HiC file (with error handling)
    JUICER_HIC(
        MERGE_DEDUP_SPLIT.out.mapped_pairs,
        chrom_sizes_ch.first()
    )
    
    // Extract QC metrics
    QC_METRICS(
        MERGE_DEDUP_SPLIT.out.stats_file
    )
    
    // Emit final outputs
    emit:
        mapped_pairs = MERGE_DEDUP_SPLIT.out.mapped_pairs
        bam = MERGE_DEDUP_SPLIT.out.bam
        bai = MERGE_DEDUP_SPLIT.out.bai
        stats = MERGE_DEDUP_SPLIT.out.stats_file
        cool = COOLER_PROCESS.out.cool
        raw_mcool = COOLER_PROCESS.out.raw_mcool
        balanced_mcool = COOLER_PROCESS.out.balanced_mcool
        hic = JUICER_HIC.out.hic
        qc_metrics = QC_METRICS.out.qc_file
}

/*
========================================================================================
    WORKFLOW COMPLETION
========================================================================================
*/

workflow.onComplete {
    log.info """
    ========================================================================================
                            Pipeline Execution Summary
    ========================================================================================
    Completed at    : ${workflow.complete}
    Duration        : ${workflow.duration}
    Success         : ${workflow.success}
    Work directory  : ${workflow.workDir}
    Exit status     : ${workflow.exitStatus}
    Error message   : ${workflow.errorMessage ?: 'None'}
    ========================================================================================
    """.stripIndent()
    
    if (workflow.success) {
        log.info "🎉 Pipeline completed successfully!"
        log.info "📁 Results are in: ${params.output_dir}"
    } else {
        log.error "💥 Pipeline failed!"
        log.error "Check the error message above and workflow logs for details."
    }
}

workflow.onError {
    log.error "Pipeline execution stopped with the following message: ${workflow.errorMessage}"
}
