// preprocess-hic: restriction-enzyme Hi-C FASTQ -> lossless BAM, MAPQ >= 1 pairs, dedup and
// restriction-fragment QC, balanced .mcool. 4DN/Open2C recipe (bwa mem -SP5M + pairtools).
// See PREPROCESS-HIC.md. Submit with preprocess-hic-submit.sh.
nextflow.enable.dsl = 2

for (p in ['sample_id', 'bwa_index', 'chrom_sizes', 'fq_input', 'outdir', 'tmpdir']) {
    if (!params[p]) error "--${p} is required"
}
if (!params.frag_bed && !params.genome_fasta) error "Set --frag_bed or --genome_fasta (for the in-silico digest)"
// An unmatched glob gives an empty channel, and Nextflow would then report success with no outputs.
if (!file(params.fq_input)) error "No FASTQs match --fq_input ${params.fq_input}"

def mcool_name = "${params.sample_id}.mapq_${params.mcool_min_mapq}.mcool"

// Split each FASTQ pair into chunks of chunk_reads reads (as microc-preprocess 1.0.3+).
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
    stem() { local b=\${1%.gz}; b=\${b%.fastq}; echo \${b%.fq}; }
    mkdir chunks/r1 chunks/r2
    mv chunks/\$(stem ${fastq_r1}).part_* chunks/r1/
    mv chunks/\$(stem ${fastq_r2}).part_* chunks/r2/
    """
}

// Restriction fragments of the genome (BED) for pairtools restrict, unless --frag_bed is given.
process digest_genome {
    input:
    path fasta
    path chrom_sizes
    output:
    path "restriction_fragments.bed"

    """
    #!/bin/bash
    set -eo pipefail
    digest_genome.py --fasta ${fasta} --chrom-sizes ${chrom_sizes} --enzymes ${params.enzymes} -o restriction_fragments.bed
    """
}

// Align a chunk, parse to pairs, and split into (a) a coordinate-sorted BAM with every alignment of every
// read (pairtools parse keeps all SAM records of a read pair, whatever its pair type) and (b) light pairs
// (no SAM columns), sorted, one bgzipped file per chrom1 for the per-chromosome dedup of the merge job.
process hic_align {
    input:
    tuple val(id), path(fastq_r1), path(fastq_r2)
    path reference_index
    path chrom_sizes
    output:
    path "${fastq_r1}.blk_*.pairs.gz", emit: pairs
    path "${fastq_r1}.blocks.tsv", emit: blocks
    path "${fastq_r1}.bam", emit: bam

    """
    #!/bin/bash
    set -eo pipefail
    trap '[ -n "${params.failed_nodes_file}" ] && hostname >> ${params.failed_nodes_file}' ERR

    export TMPDIR=${params.tmpdir}
    if [[ "${reference_index}" == *.tar.gz ]] || [[ "${reference_index}" == *.tgz ]]; then
        mkdir genome_index
        tar zxf ${reference_index} -C genome_index
        BWT=\$(find genome_index -name '*.bwt' | head -1)
        GENOME_INDEX_FA=\$(dirname "\${BWT}")/\$(basename "\${BWT}" .bwt)
    else
        # The staged symlink's real path, so that bwa finds the index files next to it.
        GENOME_INDEX_FA=\$(realpath ${reference_index})
    fi
    PARSE_NPROC=\$(( ${task.cpus} / 2 ))
    # -SP5M: no mate rescue or pairing, 5'-most alignment primary, split hits marked secondary.
    # bwa's default score threshold (-T 30). bwa-mem2 gives the same alignments as bwa mem.
    ${params.aligner} mem -SP5M -t${task.cpus} \$GENOME_INDEX_FA ${fastq_r1} ${fastq_r2} | \
    pairtools parse --min-mapq ${params.min_mapq} --walks-policy 5unique \
        --add-columns ${params.add_columns} \
        --nproc-in \${PARSE_NPROC} --nproc-out \${PARSE_NPROC} --chroms-path ${chrom_sizes} | \
    pairtools split --output-pairs unsorted.pairs.lz4 --output-sam - | \
    samtools sort -@2 -m 1G -T srt -o ${fastq_r1}.bam
    pairtools sort --nproc ${task.cpus} --memory 2G --tmpdir . unsorted.pairs.lz4 | \
    split_pairs_by_chrom1.py ${fastq_r1} 2
    rm unsorted.pairs.lz4
    """
}

// Per sample: merge the chunk BAMs (all runs) into the lossless BAM; merge and deduplicate the pairs per
// chrom1 (duplicates always share chrom1; runs and chunks are merged first), annotate restriction
// fragments, concatenate in sort order into a bgzipped, pairix-indexed .pairs.gz; QC summary.
process hic_merge {
    publishDir params.bamdir ?: params.outdir, pattern: "*.bam*", mode: 'copy', overwrite: false
    publishDir params.pairsdir ?: params.outdir, pattern: "*.pairs.gz*", mode: 'copy', overwrite: false
    publishDir params.statsdir ?: params.outdir, pattern: "*.{stats.txt,restrict_stats.tsv,flagstat.tsv,hic_qc.json}", mode: 'copy', overwrite: false
    publishDir params.outdir, pattern: "*.tool_versions.txt", mode: 'copy', overwrite: false
    input:
    path chunk_pairs
    path chunk_blocks
    path chunk_bams
    path frag_bed
    output:
    path "${params.sample_id}.pairs.gz", emit: pairs
    path "${params.sample_id}.pairs.gz.px2"
    path "${params.sample_id}.bam"
    path "${params.sample_id}.bam.bai"
    path "${params.sample_id}.dedup.stats.txt"
    path "${params.sample_id}.restrict_stats.tsv"
    path "${params.sample_id}.flagstat.tsv"
    path "${params.sample_id}.hic_qc.json"
    path "${params.sample_id}.tool_versions.txt"

    """
    #!/bin/bash
    set -eo pipefail
    trap '[ -n "${params.failed_nodes_file}" ] && hostname >> ${params.failed_nodes_file}' ERR
    SID=${params.sample_id}
    export TMPDIR=${params.tmpdir}

    # Lossless BAM: every alignment of every read, coordinate-sorted (chunk BAMs are sorted already).
    samtools merge -@4 -c -p -o \$SID.bam ${chunk_bams} &
    BAM_MERGE=\$!

    LC_ALL=C sort -t \$'\\t' -k2,2 -u ${chunk_blocks} > blocks.tsv
    mkdir blk
    cut -f1 blocks.tsv | xargs -P ${task.cpus} -I{} bash -c 'set -eo pipefail
        pairtools merge --tmpdir ${params.tmpdir} --nproc 1 *.blk_{}.pairs.gz | \
        pairtools dedup --backend ${params.dedup_backend} --output-stats blk/{}.stats | \
        pairtools restrict -f ${frag_bed} -o blk/{}.pairs
        hic_frag_stats.py blk/{}.pairs > blk/{}.frag'
    { grep '^#' blk/\$(head -1 blocks.tsv | cut -f1).pairs
      cut -f1 blocks.tsv | while read b; do grep -v '^#' blk/\$b.pairs || [ \$? = 1 ]; done
    } | bgzip -@${task.cpus} > \$SID.pairs.gz
    pairix -p pairs \$SID.pairs.gz
    pairtools stats --merge blk/*.stats -o \$SID.dedup.stats.txt
    hic_frag_stats.py --merge blk/*.frag > \$SID.restrict_stats.tsv
    rm -r blk

    wait \$BAM_MERGE
    samtools index -@${task.cpus} \$SID.bam
    samtools flagstat -@${task.cpus} -O tsv \$SID.bam > \$SID.flagstat.tsv
    hic_qc_summary.py --sample-id \$SID --version ${workflow.manifest.version} \
        --dedup-stats \$SID.dedup.stats.txt --restrict-stats \$SID.restrict_stats.tsv \
        --flagstat \$SID.flagstat.tsv -o \$SID.hic_qc.json
    { pairtools --version; { bwa 2>&1 || true; } | grep Version; echo "bwa-mem2 \$(bwa-mem2 version 2>/dev/null | tail -1)";
      samtools --version | head -1; seqkit version; echo "cooler \$(cooler --version | awk '{print \$NF}')";
      echo "pairix \$(pairix --help 2>&1 | grep -m1 -i version | awk '{print \$NF}')"; } > \$SID.tool_versions.txt
    """
}

// Contact matrix: pairs with MAPQ >= mcool_min_mapq on both sides, binned at the first resolution,
// coarsened to the others and ICE-balanced (cooler defaults) at every resolution.
process hic_mcool {
    publishDir params.coolsdir ?: params.outdir, pattern: "*.mcool", mode: 'copy', overwrite: false
    input:
    path pairs
    path chrom_sizes
    output:
    path "${params.sample_id}.mapq_${params.mcool_min_mapq}.mcool"

    """
    #!/bin/bash
    set -eo pipefail
    trap '[ -n "${params.failed_nodes_file}" ] && hostname >> ${params.failed_nodes_file}' ERR
    export TMPDIR=${params.tmpdir}
    BASE_RES=\$(echo ${params.mcool_resolutions} | cut -d, -f1)
    MCOOL=${params.sample_id}.mapq_${params.mcool_min_mapq}.mcool
    # Recorded in the cooler's metadata (read by preprocess-hic-pool.sh to refuse mixing MAPQ thresholds).
    printf '%s\\n' '{"pipeline": "preprocess-hic", "pipeline_version": "${workflow.manifest.version}", "sample_id": "${params.sample_id}", "min_mapq": ${params.mcool_min_mapq}, "pairs_filter": "(mapq1 >= ${params.mcool_min_mapq}) and (mapq2 >= ${params.mcool_min_mapq})"}' > meta.json
    bgzip -dc -@2 ${pairs} | \
    pairtools select '(mapq1 >= ${params.mcool_min_mapq}) and (mapq2 >= ${params.mcool_min_mapq})' | \
    cooler cload pairs ${params.assembly ? "--assembly ${params.assembly}" : ""} --metadata meta.json \
        -c1 2 -p1 3 -c2 4 -p2 5 ${chrom_sizes}:\$BASE_RES - base.cool
    cooler zoomify -p ${task.cpus} --balance -r ${params.mcool_resolutions} -o \$MCOOL base.cool
    rm base.cool
    """
}

workflow {
    fastq_pairs = Channel.fromFilePairs(params.fq_input, flat: true, checkIfExists: true)
    split_fastq(fastq_pairs)
    // One (id, chunk_r1, chunk_r2) tuple per chunk; sorting by name (part_NNN) keeps R1 and R2 in step.
    fastq_splits = split_fastq.out.flatMap { id, r1s, r2s ->
        r1s = [r1s].flatten().sort { it.name }
        r2s = [r2s].flatten().sort { it.name }
        if (r1s.size() != r2s.size()) error "split_fastq: ${r1s.size()} R1 but ${r2s.size()} R2 chunks for ${id}"
        [r1s, r2s].transpose().collect { r1, r2 -> [id, r1, r2] }
    }
    frags = params.frag_bed ? Channel.value(file(params.frag_bed, checkIfExists: true))
                            : digest_genome(file(params.genome_fasta, checkIfExists: true), file(params.chrom_sizes))
    hic_align(fastq_splits, file(params.bwa_index), file(params.chrom_sizes))
    hic_merge(hic_align.out.pairs.collect(), hic_align.out.blocks.collect(), hic_align.out.bam.collect(), frags)
    hic_mcool(hic_merge.out.pairs, file(params.chrom_sizes))
}

// Provenance next to the outputs: pipeline version and commit, parameters, inputs, tool versions, run metadata.
workflow.onComplete {
    def outdir = new File(params.outdir.toString())
    def sid = params.sample_id
    def vfile = new File(outdir, "${sid}.tool_versions.txt")
    def fq = file(params.fq_input)
    def fastqs = (fq instanceof List ? fq : [fq]).collect { f -> [path: f.toString(), size: f.size(), mtime: new Date(f.lastModified()).format("yyyy-MM-dd'T'HH:mm:ssXXX")] }
    def statsdir = (params.statsdir ?: params.outdir).toString()
    def qcfile = new File("${statsdir}/${sid}.hic_qc.json")
    def prov = [
        pipeline: 'preprocess-hic',
        pipeline_version: workflow.manifest.version,
        git_commit: params.git_commit,
        success: workflow.success,
        started: workflow.start.toString(),
        completed: workflow.complete.toString(),
        duration: workflow.duration.toString(),
        nextflow_version: workflow.nextflow.version.toString(),
        session_id: workflow.sessionId.toString(),
        command_line: workflow.commandLine,
        slurm_master_job: System.getenv('SLURM_JOB_ID'),
        master_node: InetAddress.getLocalHost().getHostName(),
        params: [aligner: params.aligner, bwa_mem_options: '-SP5M', min_mapq: params.min_mapq, walks_policy: '5unique',
                 add_columns: params.add_columns, dedup_backend: params.dedup_backend, dedup_by: 'chrom1',
                 chunk_reads: params.chunk_reads, enzymes: params.frag_bed ? null : params.enzymes,
                 frag_bed: params.frag_bed?.toString(), genome_fasta: params.genome_fasta?.toString(),
                 mcool_min_mapq: params.mcool_min_mapq, mcool_resolutions: params.mcool_resolutions,
                 assembly: params.assembly, bwa_index: params.bwa_index.toString(), chrom_sizes: params.chrom_sizes.toString()],
        inputs: fastqs,
        outputs: [bam: "${(params.bamdir ?: params.outdir)}/${sid}.bam",
                  pairs: "${(params.pairsdir ?: params.outdir)}/${sid}.pairs.gz",
                  dedup_stats: "${statsdir}/${sid}.dedup.stats.txt",
                  restrict_stats: "${statsdir}/${sid}.restrict_stats.tsv",
                  qc: "${statsdir}/${sid}.hic_qc.json",
                  mcool: "${(params.coolsdir ?: params.outdir)}/${mcool_name}"],
        stats_origin: 'original (written by this run)',
        tool_versions: vfile.exists() ? vfile.readLines() : [],
        qc: qcfile.exists() ? new groovy.json.JsonSlurper().parse(qcfile) : null
    ]
    outdir.mkdirs()
    new File(outdir, "${sid}.preprocess-hic.provenance.json").text = groovy.json.JsonOutput.prettyPrint(groovy.json.JsonOutput.toJson(prov)) + "\n"
}
