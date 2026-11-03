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
    .set {fastq_pairs}

// Split each FASTQ pair into chunks of chunk_reads reads in its own Slurm job. (Up to 1.0.2 Nextflow's
// splitFastq did this in the master JVM, and the page cache of the chunk writes filled the master's
// memory limit.) Chunks are gzipped at level 1 to cut work-dir I/O; bwa reads them directly.
process split_fastq {
    input:
    tuple val(id), path(fastq_r1), path(fastq_r2)
    output:
    tuple val(id), path("chunks/r1/*"), path("chunks/r2/*")

    """
    #!/bin/bash
    set -eo pipefail
    trap '[ -n "${params.failed_nodes_file}" ] && hostname >> ${params.failed_nodes_file}' ERR

    seqkit split2 -1 ${fastq_r1} -2 ${fastq_r2} -s ${params.chunk_reads} -j ${task.cpus} \
        --compress-level 1 -O chunks
    # seqkit names chunks <fastq name without .fastq/.fq[.gz]>.part_NNN.<ext>
    stem() { local b=\${1%.gz}; b=\${b%.fastq}; echo \${b%.fq}; }
    mkdir chunks/r1 chunks/r2
    mv chunks/\$(stem ${fastq_r1}).part_* chunks/r1/
    mv chunks/\$(stem ${fastq_r2}).part_* chunks/r2/
    """
}

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
    path "${fastq_r1}.blk_*.pairs.gz", emit: pairs
    path "${fastq_r1}.blocks.tsv", emit: blocks
    path "${fastq_r1}.bam", emit: bam

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
    # bwa-mem2 (default from 1.1.0) gives the same alignments as bwa mem ~2x faster; it needs the
    # bwa-mem2 index files (.0123, .bwt.2bit.64) next to the FASTA. params.aligner = 'bwa' reproduces <= 1.0.3.
    ${params.aligner} mem -5SP -T0 -t${bwa_cores} \$GENOME_INDEX_FA ${fastq_r1} ${fastq_r2} | \
    pairtools parse --min-mapq ${mapq} --walks-policy 5unique \
    --max-inter-align-gap 30 --add-columns pos5,pos3,dist_to_5,dist_to_3,read_len \
    --nproc-in \${PARSE_NPROC} --nproc-out \${PARSE_NPROC} --chroms-path ${chrom_sizes} | \
    pairtools split --output-pairs unsorted.pairs.lz4 --output-sam - | \
    samtools view -u -e '[Yt]=="UU" || [Yt]=="UR" || [Yt]=="RU"' - | \
    samtools sort -@2 -m 1G -T srt -o ${fastq_r1}.bam
    # From 1.1.0 the pairs and the alignments travel separately: light pairs (no SAM columns) to the
    # merge job's dedup, and a coordinate-sorted BAM of the chunk's alignments, kept only for pair types
    # that dedup can output (UU, UR, RU), to be merged and filtered to the deduplicated read IDs.
    # The sorted pairs are written as one file per chrom1 (bin/split_pairs_by_chrom1.py), so that the
    # merge job can merge and deduplicate the chromosomes in parallel. (lz4 and bgzip, not .gz paths:
    # for those pairtools uses pbgzip, whose threads spin on locks.)
    pairtools sort --nproc ${bwa_cores} --memory 2G --tmpdir . unsorted.pairs.lz4 | \
    split_pairs_by_chrom1.py ${fastq_r1} 2
    rm unsorted.pairs.lz4
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
    path chunk_pairs
    path chunk_blocks
    path chunk_bams
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
    # BAM: merge the chunk BAMs while the pairs are merged and deduplicated (dedup uses ~1.5 cores).
    samtools merge -@4 -c -p -l 1 -o all.bam ${chunk_bams} &
    BAM_MERGE=\$!
    # Pairs: duplicates can only share chrom1, so each chrom1 block (one file per chunk) is merged and
    # deduplicated on its own, ${task.cpus} at a time, and the outputs are concatenated in pairtools sort order
    # (chrom1 in C collation). Pairs without SAM columns: dedup is ~4x faster than on pairsam.
    # scipy (default) is transitive; its chunked clustering can see a cluster split at a different place than
    # in one whole-sample run, so a few pairs per 10 million can differ from <= 1.0.3 (UT test: 34 of 57.6M).
    LC_ALL=C sort -t \$'\t' -k2,2 -u ${chunk_blocks} > blocks.tsv
    mkdir blk
    cut -f1 blocks.tsv | xargs -P ${task.cpus} -I{} bash -c 'set -eo pipefail
        pairtools merge --tmpdir ${tmpdir} --nproc 1 *.blk_{}.pairs.gz | \
        pairtools dedup --backend ${params.dedup_backend} --mark-dups --output-stats blk/{}.stats -o blk/{}.pairs --output-dups blk/{}.dups'
    grep '^#' blk/\$(head -1 blocks.tsv | cut -f1).pairs > ${sample_id}.mapped.pairs
    cut -f1 blocks.tsv | while read b; do grep -v '^#' blk/\$b.pairs || [ \$? = 1 ]; done >> ${sample_id}.mapped.pairs
    pairtools stats --merge blk/*.stats -o ${sample_id}.stats.txt
    wait \$BAM_MERGE
    # Keep the alignments of the deduplicated pairs only: the same records as the <= 1.0.3 BAM. all.bam holds
    # UU/UR/RU pairs, each either kept or a duplicate, so keep the kept read IDs or drop the duplicates' IDs,
    # whichever list is shorter: samtools holds it in memory (~75 bytes per Illumina read name).
    grep -v '^#' ${sample_id}.mapped.pairs | cut -f1 > nodup_ids.txt
    { grep -hv '^#' blk/*.dups || [ \$? = 1 ]; } | cut -f1 > dup_ids.txt
    if [ \$(wc -l < dup_ids.txt) -lt \$(wc -l < nodup_ids.txt) ]; then IDS=^dup_ids.txt; else IDS=nodup_ids.txt; fi
    samtools view -@8 -N \$IDS -b -o ${sample_id}.bam all.bam
    rm -r all.bam nodup_ids.txt dup_ids.txt blk
    samtools index -@8 ${sample_id}.bam
    { pairtools --version; { bwa 2>&1 || true; } | grep Version; echo "bwa-mem2 \$(bwa-mem2 version 2>/dev/null | tail -1)"; samtools --version | head -1; seqkit version; } > ${sample_id}.tool_versions.txt
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
        split_fastq(fastq_pairs)
        // One (id, chunk_r1, chunk_r2) tuple per chunk; sorting by name (part_NNN) keeps R1 and R2 in step.
        fastq_splits = split_fastq.out.flatMap { id, r1s, r2s ->
            r1s = [r1s].flatten().sort { it.name }
            r2s = [r2s].flatten().sort { it.name }
            if (r1s.size() != r2s.size()) error "split_fastq: ${r1s.size()} R1 but ${r2s.size()} R2 chunks for ${id}"
            [r1s, r2s].transpose().collect { r1, r2 -> [id, r1, r2] }
        }
        microc_align(fastq_splits,params.bwa_index, params.bwa_prefix, params.chrom_sizes, 6, "16GB", "500", params.mapq, 0, params.tmpdir)
        mergepairs(params.tmpdir, params.sample_id, microc_align.out.pairs.collect(), microc_align.out.blocks.collect(), microc_align.out.bam.collect())
        printf("Completed")
}

// Provenance next to the outputs: pipeline version and commit, parameters, inputs, tool versions, run metadata.
// With params.dataset_dir set, also write a manifest row for the FOCI-3D batch pipeline (dataset_dir, pairs, sample_id, assay, stats).
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
        params: [aligner: params.aligner, dedup_backend: params.dedup_backend, dedup_by: 'chrom1', mapq: params.mapq, walks_policy: '5unique', max_inter_align_gap: 30,
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
