process EXTRACT_BWA_INDEX {
    tag "${chunk_id}"
    label 'process_low'
    
    publishDir "${params.output_dir}/align/${chunk_id}/genome_index", mode: 'copy'
    
    conda (params.enable_conda ? "bioconda::bwa=0.7.17" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bwa:0.7.17--hed695b0_7' :
        'quay.io/biocontainers/bwa:0.7.17--hed695b0_7' }"
    
    input:
    val chunk_id
    path bwa_index_tar
    
    output:
    tuple val(chunk_id), path("genome_index"), emit: genome_index
    path "versions.yml", emit: versions
    
    when:
    task.ext.when == null || task.ext.when
    
    script:
    def args = task.ext.args ?: ''
    """
    # Create genome index directory
    mkdir -p genome_index
    
    # Extract BWA index (compatible with BusyBox tar)
    if tar --help 2>&1 | grep -q "GNU tar"; then
        # GNU tar supports -z flag
        tar zxvf ${bwa_index_tar} -C genome_index
    else
        # BusyBox tar doesn't support -z, use gzip -dc + tar
        gzip -dc ${bwa_index_tar} | tar xvf - -C genome_index
    fi
    
    # Find the reference FASTA file (should have .bwt extension companion)
    bwt_file=\$(find genome_index -name "*.bwt" | head -1)
    if [ -z "\$bwt_file" ]; then
        echo "Error: No BWA index files found in ${bwa_index_tar}"
        exit 1
    fi
    
    # Log the extracted files
    echo "Extracted BWA index files:"
    ls -la genome_index/
    
    # Create versions file
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwa: \$(echo \$(bwa 2>&1) | sed 's/^.*Version: //; s/Contact:.*\$//')
        tar: \$(tar --version 2>/dev/null | head -1 | sed 's/tar (GNU tar) //' || echo "BusyBox tar")
    END_VERSIONS
    """
    
    stub:
    """
    mkdir -p genome_index
    touch genome_index/test.fa.gz
    touch genome_index/test.fa.gz.amb
    touch genome_index/test.fa.gz.ann
    touch genome_index/test.fa.gz.bwt
    touch genome_index/test.fa.gz.pac
    touch genome_index/test.fa.gz.sa
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwa: 0.7.17
        tar: 1.34
    END_VERSIONS
    """
}
