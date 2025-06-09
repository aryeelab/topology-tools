process MERGE_DEDUP_SPLIT {
    tag "${sample_id}"
    label 'process_high'
    
    publishDir "${params.output_dir}/merge", mode: 'copy'
    
    conda (params.enable_conda ? "bioconda::pairtools=1.0.2 bioconda::samtools=1.15" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-pairtools-samtools:latest' :
        'quay.io/biocontainers/mulled-v2-pairtools-samtools:latest' }"
    
    input:
    tuple val(sample_id), path(pairsam_files)
    
    output:
    tuple val(sample_id), path("*.mapped.pairs"), emit: mapped_pairs
    tuple val(sample_id), path("*.stats.txt"), emit: stats_file
    tuple val(sample_id), path("*.bam"), emit: bam
    tuple val(sample_id), path("*.bam.bai"), emit: bai
    path "versions.yml", emit: versions
    
    when:
    task.ext.when == null || task.ext.when
    
    script:
    def args = task.ext.args ?: ''
    def mapped_pairs = "${sample_id}.mapped.pairs"
    def stats_file = "${sample_id}.stats.txt"
    def bam_file = "${sample_id}.bam"
    def merge_cores = task.cpus > 12 ? 12 : task.cpus
    def dedup_cores_in = task.cpus > 2 ? 2 : 1
    def dedup_cores_out = task.cpus > 8 ? 8 : task.cpus
    def samtools_cores = task.cpus > 6 ? 6 : task.cpus
    """
    echo "Processing ${pairsam_files.size()} pairsam files for sample ${sample_id}"
    echo "Input files:"
    for file in ${pairsam_files.join(' ')}; do
        echo "  \$file (\$(du -h \$file | cut -f1))"
    done
    
    echo "Starting merge, dedup, and split process..."
    
    # Merge, deduplicate, and split pairs
    pairtools merge \\
        --nproc ${merge_cores} \\
        ${pairsam_files.join(' ')} | \\
    pairtools dedup \\
        --nproc-in ${dedup_cores_in} \\
        --nproc-out ${dedup_cores_out} \\
        --mark-dups \\
        --output-stats ${stats_file} | \\
    pairtools split \\
        --nproc-in ${dedup_cores_in} \\
        --nproc-out ${dedup_cores_out} \\
        --output-pairs ${mapped_pairs} \\
        --output-sam - | \\
    samtools view -bS -@${samtools_cores} | \\
    samtools sort -@${samtools_cores} -o ${bam_file}
    
    # Index the BAM file
    samtools index ${bam_file}
    
    # Verify output files were created
    for file in ${mapped_pairs} ${stats_file} ${bam_file} ${bam_file}.bai; do
        if [ ! -f "\$file" ]; then
            echo "Error: Output file \$file was not created"
            exit 1
        fi
        echo "Created \$file (\$(du -h \$file | cut -f1))"
    done
    
    # Show some statistics
    echo "Statistics summary:"
    echo "  Mapped pairs: \$(wc -l < ${mapped_pairs} | sed 's/^[[:space:]]*//')"
    echo "  BAM file size: \$(du -h ${bam_file} | cut -f1)"
    
    # Create versions file
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pairtools: \$(pairtools --version 2>&1 | sed 's/pairtools, version //')
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
    
    stub:
    def mapped_pairs = "${sample_id}.mapped.pairs"
    def stats_file = "${sample_id}.stats.txt"
    def bam_file = "${sample_id}.bam"
    """
    touch ${mapped_pairs}
    touch ${stats_file}
    touch ${bam_file}
    touch ${bam_file}.bai
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pairtools: 1.0.2
        samtools: 1.15
    END_VERSIONS
    """
}
