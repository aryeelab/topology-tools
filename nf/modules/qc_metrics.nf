process QC_METRICS {
    tag "${sample_id}"
    label 'process_low'
    
    publishDir "${params.output_dir}/qc", mode: 'copy'
    
    conda (params.enable_conda ? "conda-forge::coreutils=9.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
        'ubuntu:20.04' }"
    
    input:
    tuple val(sample_id), path(stats_file)
    
    output:
    tuple val(sample_id), path("*_qc.txt"), emit: qc_file
    tuple val(sample_id), path("*_qc_summary.json"), emit: qc_json
    path "versions.yml", emit: versions
    
    when:
    task.ext.when == null || task.ext.when
    
    script:
    def args = task.ext.args ?: ''
    def qc_file = "${sample_id}_qc.txt"
    def qc_json = "${sample_id}_qc_summary.json"
    """
    echo "Extracting QC metrics for sample ${sample_id}"
    echo "Input stats file: ${stats_file} (\$(du -h ${stats_file} | cut -f1))"
    
    # Check if stats file exists and is not empty
    if [ ! -s "${stats_file}" ]; then
        echo "Error: Stats file ${stats_file} is empty or does not exist"
        exit 1
    fi
    
    echo "Contents of stats file:"
    head -20 ${stats_file}
    echo ""
    
    # Extract QC metrics in the same order as Python pipeline
    echo "Extracting QC metrics..."
    
    # Extract total reads
    total_reads=\$(cat ${stats_file} | grep -w "total" | cut -f2 || echo "0")
    echo "Total reads: \$total_reads"
    
    # Extract mapped reads
    mapped_reads=\$(cat ${stats_file} | grep -w "total_mapped" | cut -f2 || echo "0")
    echo "Mapped reads: \$mapped_reads"
    
    # Extract deduplicated reads
    nodups_reads=\$(cat ${stats_file} | grep -w "total_nodups" | cut -f2 || echo "0")
    echo "No-dups reads: \$nodups_reads"
    
    # Extract cis reads (1kb+)
    cis_1kb_reads=\$(cat ${stats_file} | grep -w "cis_1kb+" | cut -f2 || echo "0")
    echo "Cis 1kb+ reads: \$cis_1kb_reads"
    
    # Extract cis reads (10kb+)
    cis_10kb_reads=\$(cat ${stats_file} | grep -w "cis_10kb+" | cut -f2 || echo "0")
    echo "Cis 10kb+ reads: \$cis_10kb_reads"
    
    # Write QC metrics to text file (same format as Python pipeline)
    cat > ${qc_file} <<EOF
\$total_reads
\$mapped_reads
\$nodups_reads
\$cis_1kb_reads
\$cis_10kb_reads
EOF
    
    echo "Created QC file: ${qc_file}"
    cat ${qc_file}
    
    # Calculate percentages for JSON summary
    if [ "\$total_reads" -gt 0 ]; then
        mapped_pct=\$(echo "scale=2; \$mapped_reads * 100 / \$total_reads" | bc -l 2>/dev/null || echo "0")
        nodups_pct=\$(echo "scale=2; \$nodups_reads * 100 / \$total_reads" | bc -l 2>/dev/null || echo "0")
        cis_1kb_pct=\$(echo "scale=2; \$cis_1kb_reads * 100 / \$total_reads" | bc -l 2>/dev/null || echo "0")
        cis_10kb_pct=\$(echo "scale=2; \$cis_10kb_reads * 100 / \$total_reads" | bc -l 2>/dev/null || echo "0")
    else
        mapped_pct="0"
        nodups_pct="0"
        cis_1kb_pct="0"
        cis_10kb_pct="0"
    fi
    
    # Create JSON summary
    cat > ${qc_json} <<EOF
{
    "sample_id": "${sample_id}",
    "qc_metrics": {
        "reads_total": "\$total_reads",
        "reads_mapped": "\$mapped_reads",
        "reads_nodups": "\$nodups_reads",
        "reads_cis_1kb": "\$cis_1kb_reads",
        "reads_cis_10kb": "\$cis_10kb_reads"
    },
    "percentages": {
        "mapped_percent": "\$mapped_pct",
        "nodups_percent": "\$nodups_pct",
        "cis_1kb_percent": "\$cis_1kb_pct",
        "cis_10kb_percent": "\$cis_10kb_pct"
    },
    "quality_flags": {
        "sufficient_reads": \$([ "\$total_reads" -gt 1000000 ] && echo "true" || echo "false"),
        "good_mapping": \$([ "\$mapped_pct" != "0" ] && [ \$(echo "\$mapped_pct > 50" | bc -l 2>/dev/null || echo 0) -eq 1 ] && echo "true" || echo "false"),
        "low_duplication": \$([ "\$nodups_pct" != "0" ] && [ \$(echo "\$nodups_pct > 70" | bc -l 2>/dev/null || echo 0) -eq 1 ] && echo "true" || echo "false")
    }
}
EOF
    
    echo "Created QC JSON summary: ${qc_json}"
    cat ${qc_json}
    
    echo "QC metrics extraction completed successfully"
    
    # Create versions file
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        grep: \$(grep --version | head -1 | sed 's/grep (GNU grep) //')
        cut: \$(cut --version | head -1 | sed 's/cut (GNU coreutils) //')
        bc: \$(bc --version | head -1 | sed 's/bc //')
    END_VERSIONS
    """
    
    stub:
    def qc_file = "${sample_id}_qc.txt"
    def qc_json = "${sample_id}_qc_summary.json"
    """
    echo "100000" > ${qc_file}
    echo "80000" >> ${qc_file}
    echo "75000" >> ${qc_file}
    echo "50000" >> ${qc_file}
    echo "30000" >> ${qc_file}
    
    cat > ${qc_json} <<EOF
{
    "sample_id": "${sample_id}",
    "qc_metrics": {
        "reads_total": "100000",
        "reads_mapped": "80000",
        "reads_nodups": "75000",
        "reads_cis_1kb": "50000",
        "reads_cis_10kb": "30000"
    }
}
EOF
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        grep: 3.7
        cut: 8.32
        bc: 1.07.1
    END_VERSIONS
    """
}
