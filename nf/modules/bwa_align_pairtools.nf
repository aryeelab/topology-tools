process BWA_ALIGN_PAIRTOOLS {
    tag "${chunk_id}"
    label 'process_high'
    
    publishDir "${params.output_dir}/align/${chunk_id}", mode: 'copy'
    
    conda (params.enable_conda ? "bioconda::bwa=0.7.17 bioconda::pairtools=1.0.2 bioconda::samtools=1.15" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-bwa-pairtools-samtools:latest' :
        'quay.io/biocontainers/mulled-v2-bwa-pairtools-samtools:latest' }"
    
    input:
    tuple val(chunk_id), path(fastq_r1), path(fastq_r2)
    tuple val(chunk_id_index), path(genome_index_dir)
    path chrom_sizes
    
    output:
    tuple val(chunk_id), path("*.pairsam.gz"), emit: pairsam
    path "versions.yml", emit: versions
    
    when:
    task.ext.when == null || task.ext.when
    
    script:
    def args = task.ext.args ?: ''
    def mapq = params.mapq ?: 20
    def bwa_cores = params.bwa_cores ?: 2
    def output_pairsam = "${fastq_r1.baseName}.pairsam.gz"
    """
    # Find the reference genome file (follow symlinks with -L)
    bwt_file=\$(find -L ${genome_index_dir} -name "*.bwt" | head -1)
    if [ -z "\$bwt_file" ]; then
        echo "Error: No BWA index files found in ${genome_index_dir}"
        exit 1
    fi
    
    # Get the reference genome prefix (remove .bwt extension)
    genome_index_fa=\${bwt_file%.bwt}
    
    echo "Using BWA index: \$genome_index_fa"
    echo "Processing FASTQ files: ${fastq_r1}, ${fastq_r2}"
    echo "Output file: ${output_pairsam}"
    
    # Run BWA alignment piped to pairtools
    bwa mem -5SP -T0 -t${bwa_cores} \$genome_index_fa ${fastq_r1} ${fastq_r2} | \\
    pairtools parse \\
        --min-mapq ${mapq} \\
        --walks-policy 5unique \\
        --max-inter-align-gap 30 \\
        --add-columns pos5,pos3,dist_to_5,dist_to_3,read_len \\
        --nproc-in ${bwa_cores} \\
        --nproc-out ${bwa_cores} \\
        --chroms-path ${chrom_sizes} | \\
    pairtools sort \\
        --nproc ${bwa_cores} \\
        -o ${output_pairsam}
    
    # Verify output file was created
    if [ ! -f "${output_pairsam}" ]; then
        echo "Error: Output file ${output_pairsam} was not created"
        exit 1
    fi
    
    echo "Successfully created ${output_pairsam}"
    echo "File size: \$(du -h ${output_pairsam})"
    
    # Create versions file
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwa: \$(echo \$(bwa 2>&1) | sed 's/^.*Version: //; s/Contact:.*\$//')
        pairtools: \$(pairtools --version 2>&1 | sed 's/pairtools, version //')
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
    
    stub:
    def output_pairsam = "${fastq_r1.baseName}.pairsam.gz"
    """
    touch ${output_pairsam}
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwa: 0.7.17
        pairtools: 1.0.2
        samtools: 1.15
    END_VERSIONS
    """
}
