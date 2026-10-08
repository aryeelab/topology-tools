process juicer_hic {
    publishDir params.outdir
    input:
    path outdir
    val sample_id
    path chrom_sizes
    path mapped_pairs
    val cores
    output:
    path "${sample_id}.hic"

    """
    #!/bin/bash
    set -eo pipefail
    trap '[ -n "${params.failed_nodes_file}" ] && hostname >> ${params.failed_nodes_file}' ERR

    java -Xmx120g  -Djava.awt.headless=true -jar /usr/local/bin/juicer_tools_1.22.01.jar pre \
        --threads ${cores} \
        ${outdir}/${mapped_pairs} \
        ${sample_id}.hic \
        ${chrom_sizes}
    """
}

process cooler {
    publishDir params.outdir
    input:
    path outdir
    val sample_id
    path chrom_sizes
    path mapped_pairs
    val resolution
    output:
    path "${outdir}/${sample_id}.raw.mcool"   
    path "${outdir}/${sample_id}.balanced.mcool"  
    
    """
    #!/bin/bash
    set -eo pipefail
    trap '[ -n "${params.failed_nodes_file}" ] && hostname >> ${params.failed_nodes_file}' ERR

    cooler cload pairs -c1 2 -p1 3 -c2 4 -p2 5 ${chrom_sizes}:250 ${outdir}/${mapped_pairs} ${outdir}/${sample_id}.cool
    cooler zoomify --resolutions ${resolution} -o ${outdir}/${sample_id}.raw.mcool -p 4 ${outdir}/${sample_id}.cool
    cooler zoomify --resolutions ${resolution} -o ${outdir}/${sample_id}.balanced.mcool -p 4 --balance --balance-args '--nproc 4' ${outdir}/${sample_id}.cool
    """
}

process footprint_counts {
    publishDir params.outdir
    input:
    path outdir
    val sample_id
    path mapped_pairs
    output:
    path "${sample_id}.counts.tsv.gz"
    """
    #!/bin/bash
    set -eo pipefail
    trap '[ -n "${params.failed_nodes_file}" ] && hostname >> ${params.failed_nodes_file}' ERR

    pairs_to_fragment_counts.py ${mapped_pairs} -o ${sample_id}.counts.tsv.gz
    """
}

process qc {
    publishDir params.outdir
    input:
    path outdir
    path mapped_stats
    val sample_id
    output:
    path "${sample_id}.qc.txt"
    
    """
    #!/bin/bash
    set -eo pipefail
    trap '[ -n "${params.failed_nodes_file}" ] && hostname >> ${params.failed_nodes_file}' ERR
    cat ${mapped_stats} | grep -w "total" | cut -f2 > ${sample_id}.qc.txt
    cat ${mapped_stats} | grep -w "total_mapped" | cut -f2 > ${sample_id}.qc.txt
    cat ${mapped_stats} | grep -w "total_nodups" | cut -f2 > ${sample_id}.qc.txt
    cat ${mapped_stats} | grep -w "cis_1kb+" | cut -f2 > ${sample_id}.qc.txt
    cat ${mapped_stats} | grep -w "cis_10kb+" | cut -f2 > ${sample_id}.qc.txt
    """
}

process qc_py {
    publishDir params.outdir
    input:
    path outdir
    path bam
    path mapped_pairs
    path converter_script
    val sample_id
    output:
    path "${sample_id}.qc.json"
    path "${sample_id}.fragment_lengths.tsv"
    path "${sample_id}.sample_metrics.tsv"
    path "${sample_id}.fragment_lengthsnoheader.tsv"
    """
    #!/bin/bash
    set -eo pipefail
    trap '[ -n "${params.failed_nodes_file}" ] && hostname >> ${params.failed_nodes_file}' ERR
    microc-qc --pairs ${mapped_pairs} --bam ${bam} --out ${sample_id}.qc.json
    bash ${converter_script} -f ${sample_id}.fragment_lengths.tsv -m ${sample_id}.sample_metrics.tsv ${sample_id}.qc.json
    tail -n+2 ${sample_id}.fragment_lengths.tsv > ${sample_id}.fragment_lengthsnoheader.tsv
    """
}
workflow {
    main:
        juicer_hic(params.outdir, params.sample_id, params.chrom_sizes, params.mapped_pairs, "4")
        cooler(params.outdir, params.sample_id, params.chrom_sizes, params.mapped_pairs, "500,1000,2000,5000,10000,20000")
        footprint_counts(params.outdir, params.sample_id, params.mapped_pairs)
        qc(params.outdir, params.mapped_stats, params.sample_id)
        qc_py(params.outdir, params.bam, params.mapped_pairs, "/cluster/aryeelab/mark/convert_microc_json.sh", params.sample_id)
}
