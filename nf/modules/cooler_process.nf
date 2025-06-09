process COOLER_PROCESS {
    tag "${sample_id}"
    label 'process_medium'
    
    publishDir "${params.output_dir}/cooler", mode: 'copy'
    
    conda (params.enable_conda ? "bioconda::cooler=0.9.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/cooler:0.9.1--pyh7cba7a3_0' :
        'quay.io/biocontainers/cooler:0.9.1--pyh7cba7a3_0' }"
    
    input:
    tuple val(sample_id), path(mapped_pairs)
    path chrom_sizes
    
    output:
    tuple val(sample_id), path("*.cool"), emit: cool
    tuple val(sample_id), path("*.raw.mcool"), emit: raw_mcool
    tuple val(sample_id), path("*.balanced.mcool"), emit: balanced_mcool
    path "versions.yml", emit: versions
    
    when:
    task.ext.when == null || task.ext.when
    
    script:
    def args = task.ext.args ?: ''
    def resolution = params.resolution ?: 1000
    def cool_file = "${sample_id}.cool"
    def raw_mcool = "${sample_id}.raw.mcool"
    def balanced_mcool = "${sample_id}.balanced.mcool"
    def cooler_cores = task.cpus > 4 ? 4 : task.cpus
    """
    echo "Processing cooler files for sample ${sample_id}"
    echo "Input mapped pairs: ${mapped_pairs} (\$(du -h ${mapped_pairs} | cut -f1))"
    echo "Chromosome sizes: ${chrom_sizes}"
    echo "Resolution: ${resolution} bp"
    
    # Generate .cool file
    echo "Step 1: Generating .cool file..."
    cooler cload pairs \\
        -c1 2 -p1 3 -c2 4 -p2 5 \\
        ${chrom_sizes}:${resolution} \\
        ${mapped_pairs} \\
        ${cool_file}
    
    if [ ! -f "${cool_file}" ]; then
        echo "Error: Failed to create ${cool_file}"
        exit 1
    fi
    echo "Created ${cool_file} (\$(du -h ${cool_file} | cut -f1))"
    
    # Generate raw mcool file (multi-resolution)
    echo "Step 2: Generating raw mcool file..."
    cooler zoomify \\
        --resolutions ${resolution}N \\
        -o ${raw_mcool} \\
        -p ${cooler_cores} \\
        ${cool_file}
    
    if [ ! -f "${raw_mcool}" ]; then
        echo "Error: Failed to create ${raw_mcool}"
        exit 1
    fi
    echo "Created ${raw_mcool} (\$(du -h ${raw_mcool} | cut -f1))"
    
    # Generate balanced mcool file (with matrix balancing)
    echo "Step 3: Generating balanced mcool file..."
    cooler zoomify \\
        --resolutions ${resolution}N \\
        -o ${balanced_mcool} \\
        -p ${cooler_cores} \\
        --balance \\
        --balance-args '--nproc ${cooler_cores}' \\
        ${cool_file}
    
    if [ ! -f "${balanced_mcool}" ]; then
        echo "Error: Failed to create ${balanced_mcool}"
        exit 1
    fi
    echo "Created ${balanced_mcool} (\$(du -h ${balanced_mcool} | cut -f1))"
    
    # Show summary
    echo "Cooler processing completed successfully:"
    echo "  Cool file: ${cool_file} (\$(du -h ${cool_file} | cut -f1))"
    echo "  Raw mcool: ${raw_mcool} (\$(du -h ${raw_mcool} | cut -f1))"
    echo "  Balanced mcool: ${balanced_mcool} (\$(du -h ${balanced_mcool} | cut -f1))"
    
    # Create versions file
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cooler: \$(cooler --version 2>&1 | sed 's/cooler, version //')
    END_VERSIONS
    """
    
    stub:
    def cool_file = "${sample_id}.cool"
    def raw_mcool = "${sample_id}.raw.mcool"
    def balanced_mcool = "${sample_id}.balanced.mcool"
    """
    touch ${cool_file}
    touch ${raw_mcool}
    touch ${balanced_mcool}
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cooler: 0.9.1
    END_VERSIONS
    """
}
