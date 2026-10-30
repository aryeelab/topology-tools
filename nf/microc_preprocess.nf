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

// An unmatched glob gives an empty channel, and Nextflow would then report success with no outputs.
if (!file(params.fq_input)) error "No FASTQs match --fq_input ${params.fq_input}"

Channel
    .fromFilePairs(params.fq_input, flat:true, checkIfExists:true)
    .splitFastq(by: params.chunk_reads, pe:true, file:true)
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
    if [[ "${reference_index}" == *.tar.gz ]] || [[ "${reference_index}" == *.tgz ]]; then
        echo "Using provided reference .tar.gz: ${reference_index}"
        outdir="genome_index"
        mkdir \${outdir}
        tar zxvf ${reference_index} -C \${outdir}
        BWT=\$(find "\${outdir}" -name '*.bwt')
        GENOME_INDEX_FA=\$(dirname "\${BWT}")/\$(basename "\${BWT}" .bwt)
    elif [ ! -z "${reference_index}" ]; then
        echo "Using filesystem index prefix: ${reference_index}"
        # Follow the Nextflow-staged symlink to get the real path so BWA can
        # find the companion .amb/.ann/.bwt/.pac/.sa files alongside it
        GENOME_INDEX_FA=\$(realpath ${reference_index})
    else
        echo "Using GCS reference_index_prefix: ${reference_index_prefix}"
        mkdir genome_index
        cd genome_index
        gsutil cp ${reference_index_prefix}.amb .
        gsutil cp ${reference_index_prefix}.ann .
        gsutil cp ${reference_index_prefix}.bwt .
        gsutil cp ${reference_index_prefix}.pac .
        gsutil cp ${reference_index_prefix}.sa .
        cd ..
        BWT=\$(find "genome_index" -name '*.bwt')
        GENOME_INDEX_FA=\$(dirname "\${BWT}")/\$(basename "\${BWT}" .bwt)
    fi
    echo "Using BWA index: \$GENOME_INDEX_FA"
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
    // mode 'copy': the default ('symlink') leaves outputs as links into the Nextflow work dir on scratch.
    // overwrite: false: never replace existing lab pairs/BAM/stats.
    publishDir params.bamdir ?: params.outdir, pattern: "*.bam*", mode: 'copy', overwrite: false
    publishDir params.pairsdir ?: params.outdir, pattern: "*.mapped.pairs", mode: 'copy', overwrite: false
    publishDir params.statsdir ?: params.outdir, pattern: "*.stats.txt", mode: 'copy', overwrite: false
    publishDir params.outdir, pattern: "*.tool_versions.txt", mode: 'copy', overwrite: false
    input:
    path tmpdir
    val sample_id
    path pairsams
    output:
    path "${sample_id}.mapped.pairs"
    // path microc_stats
    path "${sample_id}.bam"
    path "${sample_id}.bam.bai"
    path "${sample_id}.stats.txt"
    path "${sample_id}.tool_versions.txt"
    """
    #!/bin/bash
    set -eo pipefail
    trap '[ -n "${params.failed_nodes_file}" ] && hostname >> ${params.failed_nodes_file}' ERR

    export TMPDIR=${tmpdir}
    echo ${pairsams}
    pairtools merge --tmpdir ${tmpdir} -o merge.pairs.gz --nproc 4 ${pairsams}
    pairtools dedup merge.pairs.gz --nproc-in 4 --nproc-out 4 --mark-dups --output-stats ${sample_id}.stats.txt | \
    pairtools split --nproc-in 4 --nproc-out 4 --output-pairs ${sample_id}.mapped.pairs --output-sam -| \
    samtools view -bS -@4 | \
    samtools sort -@4 -o ${sample_id}.bam

    samtools index ${sample_id}.bam
    { pairtools --version; { bwa 2>&1 || true; } | grep Version; samtools --version | head -1; } > ${sample_id}.tool_versions.txt
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
        microc_align(fastq_splits,params.bwa_index, params.bwa_prefix, params.chrom_sizes, 6, "16GB", "500", params.mapq, 0, params.tmpdir)
        mergepairs(params.tmpdir, params.sample_id, microc_align.out.collect())
        printf("Completed")
}

// Provenance next to the outputs: pipeline version and commit, parameters, inputs, tool versions, run metadata.
// With params.dataset_dir set, also write a manifest row for foci-3d-dfci run_dataset.sh.
workflow.onComplete {
    def outdir = new File(params.outdir.toString())
    def sid = params.sample_id
    def vfile = new File(outdir, "${sid}.tool_versions.txt")
    def pairsdir = (params.pairsdir ?: params.outdir).toString()
    def statsdir = (params.statsdir ?: params.outdir).toString()
    def fq = file(params.fq_input)
    def fastqs = (fq instanceof List ? fq : [fq]).collect { f -> [path: f.toString(), size: f.size(), mtime: new Date(f.lastModified()).format("yyyy-MM-dd'T'HH:mm:ssXXX")] }
    def prov = [
        pipeline: 'microc-preprocess',
        pipeline_version: workflow.manifest.version,
        git_commit: params.git_commit,
        success: workflow.success,
        started: workflow.start.toString(),
        completed: workflow.complete.toString(),
        duration: workflow.duration.toString(),
        nextflow_version: workflow.nextflow.version.toString(),
        session_id: workflow.sessionId.toString(),
        command_line: workflow.commandLine,
        params: [mapq: params.mapq, walks_policy: '5unique', max_inter_align_gap: 30,
                 add_columns: 'pos5,pos3,dist_to_5,dist_to_3,read_len', chunk_reads: params.chunk_reads,
                 bwa_index: params.bwa_index.toString(), chrom_sizes: params.chrom_sizes.toString()],
        inputs: fastqs,
        outputs: [pairs: "${pairsdir}/${sid}.mapped.pairs", bam: "${(params.bamdir ?: params.outdir)}/${sid}.bam",
                  stats: "${statsdir}/${sid}.stats.txt"],
        tool_versions: vfile.exists() ? vfile.readLines() : []
    ]
    outdir.mkdirs()
    new File(outdir, "${sid}.preprocess.provenance.json").text = groovy.json.JsonOutput.prettyPrint(groovy.json.JsonOutput.toJson(prov)) + "\n"
    if (workflow.success && params.dataset_dir && new File("${pairsdir}/${sid}.mapped.pairs").exists()) {
        new File(outdir, "${sid}.foci_manifest.tsv").text =
            [params.dataset_dir, "${pairsdir}/${sid}.mapped.pairs", sid, params.assay, "${statsdir}/${sid}.stats.txt"].join('\t') + "\n"
    }
}
