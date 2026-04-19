// nextflow.enable.dsl=2
 
// include { sort } from './sort.nf'
// include { split } from './split.nf'
// include { merge } from './merge.nf'

// process sum_fastq_size {
//     input:
//         path fastq_r1
//         path fastq_r2
//     output:
//         val size_gb
    
//     """
//     #!/usr/bin/python

//     size_gb = round(size(fastq_r1, "GB") + size(fastq_r2, "GB"))
//     """
// }

// chmod +x ${top_monitor_script}
    // ${top_monitor_script} > ${chunk_id}.top.log &

    // # Start resource monitoring script
    // chmod +x ${resource_monitor_script}
    // ${resource_monitor_script} > ${chunk_id}.resources.log &

    // # Check if bwa index is provided as a tar.gz file ("reference_index")
    // #  or a URI prefix ("reference_index_prefix", e.g. gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta.64)
    // if [ ! -z "${reference_index}" ] && [ ! -z "${reference_index_prefix}" ]
    // then
    //     echo "ERROR: Both reference_index and reference_index_prefix provided. Please provide only one."
    //     exit
    // fi

Channel
    .fromFilePairs(params.fq_input, flat:true)
    .splitFastq(by: 30000000, pe:true, file:true)
    .set {fastq_splits}
process microc_align {
    // clusterOptions '-x node26,node20,node14'
    input:
    // val image_id
    // val sample_id
    // path fastq_r1
    // path fastq_r2
    tuple val(id), path(fastq_r1), path(fastq_r2)
    //tuple val(id), path(fastqs)
    path reference_index
    val reference_index_prefix
    path chrom_sizes
    val bwa_cores
    val mem
    val dis
    val mapq
    val preemptible
    path tmpdir
    // val chunk_id
    // path resource_monitor_script, name:
    // path top_monitor_script, name:
    output:
    path "${fastq_r1}.pairsam.gz"

    """
    #!/bin/bash
    set -eo pipefail
    trap '[ -n "${params.failed_nodes_file}" ] && hostname >> ${params.failed_nodes_file}' ERR

    export TMPDIR=${tmpdir}
    echo ${fastq_r1} ${fastq_r2}
    if [ ! -z "${reference_index}" ]
    then
        echo "Using provided reference .tar.gz: ${reference_index}"
        outdir="genome_index"
        mkdir \${outdir}
        tar zxvf ${reference_index} -C \${outdir}
    else
        echo "Using reference_index_prefix: ${reference_index_prefix}"
        mkdir genome_index
        cd genome_index
        gsutil cp ${reference_index_prefix}.amb .
        gsutil cp ${reference_index_prefix}.ann .
        gsutil cp ${reference_index_prefix}.bwt .
        gsutil cp ${reference_index_prefix}.pac .
        gsutil cp ${reference_index_prefix}.sa .
        echo "Downloaded bwa index files:"
        ls -lh
        cd ..
    fi

        # Get genome index name
    BWT=`find "\$outdir" -name '*.bwt'`
    GENOME_INDEX_FA=`dirname "\$BWT"`/`basename "\$BWT" .bwt`
    echo "Using bwa index: \$GENOME_INDEX_FA"
    PARSE_NPROC=\$(( ${bwa_cores} / 2 ))
    bwa mem -5SP -T0 -t${bwa_cores} \$GENOME_INDEX_FA ${fastq_r1} ${fastq_r2} | \
    pairtools parse --min-mapq ${mapq} --walks-policy 5unique \
    --max-inter-align-gap 30 --add-columns pos5,pos3,dist_to_5,dist_to_3,read_len \
    --nproc-in \${PARSE_NPROC} --nproc-out \${PARSE_NPROC} --chroms-path ${chrom_sizes} | \
    pairtools sort --nproc ${bwa_cores} -o ${fastq_r1}.pairsam.gz
    """
}
process mergepairs {
    // clusterOptions '-x node26,node20,node14'
    publishDir params.outdir
    input:
    path tmpdir
    val sample_id
    path pairsams
    output:
    path "${sample_id}.mapped.pairs"
    // path microc_stats
    path "${sample_id}.bam"
    path "${sample_id}.bam.bai"
    """
    #!/bin/bash
    set -eo pipefail
    trap '[ -n "${params.failed_nodes_file}" ] && hostname >> ${params.failed_nodes_file}' ERR

    export TMPDIR=${tmpdir}
    echo ${pairsams}
    pairtools merge --tmpdir ${tmpdir} -o merge.pairs.gz --nproc 4 sep=' ' ${pairsams}
    pairtools dedup merge.pairs.gz --nproc-in 4 --nproc-out 4 --mark-dups --output-stats ${sample_id}.stats.txt | \
    pairtools split --nproc-in 4 --nproc-out 4 --output-pairs ${sample_id}.mapped.pairs --output-sam -| \
    samtools view -bS -@4 | \
    samtools sort -@4 -o ${sample_id}.bam

    samtools index ${sample_id}.bam
    """
}

// workflow {
//     take:
//     image_id
//     pipeline_version
//     sample_id
//     fastq_r1
//     fastq_r2
//     reference_bwa_idx
//     reference_bwa_idx_prefix
//     chrom_sizes
//     num_reads_per_chunk
//     resource_monitor_script
//     top_monitor_script
//     main:
//     split(input_ch)
//     sort(split.out.flatten())
//     merge(sort.out.collect())
//     sum_fastq_size(fastq_r1,fastq_r2)
//     microc_align()
//     merge_pairs()
// }
workflow {
    // take:
    // fastq_r1
    // fastq_r2
    main:
        microc_align(fastq_splits,params.bwa_index, params.bwa_prefix, params.chrom_sizes, 6, "16GB", "500", "20", 0, params.tmpdir)
        mergepairs(params.tmpdir, params.sample_id, microc_align.out.collect())
        printf("Completed")
}
