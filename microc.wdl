version 1.0


workflow microc {
    String pipeline_ver = 'dev'
    String image_id = sub(pipeline_ver, "dev", "latest")
    
    meta {

        author: ''
        email: ''
        description: 'Process microC data. Inspired in pipeline from: https://micro-c.readthedocs.io/en/latest/index.html'
        organization: ''

        parameter_group: {
            reference_genome: {
                title: 'Reference genome',
                description: 'Genome specific files.',
                help: 'For now, only working with hardcoded hg38 genome. The idea is to input the reference genome with the indexes'
            },
            input_genomic_data: {
                title: 'Input genomic data',
                description: 'Genomic input files for experiment.',
                help: 'Paired end fastq microC files'
            },
            alignment: {
                title: 'Alignment',
                description: 'Parameters for alignment.',
                help: 'The alignment is done using bwa, and requires the reference genome indexes, tar and gziped. The number of cores is hardoded for now.'
            }
        }
    }

    input {
        String sample_id
        File fastq_R1
        File fastq_R2
        File? reference_bwa_idx
        String? reference_bwa_idx_prefix
        File chrom_sizes
        Int num_reads_per_chunk = 10000000
        File resource_monitor_script
        File top_monitor_script
    }

    # Split the comma-separated string of fastq file names into an array
    call split_string_into_array as fastq1 {input: str = fastq_R1}
    call split_string_into_array as fastq2 {input: str = fastq_R2}
    
    # Calculate the total fastq file size
    #scatter (fq1 in fastq1.out) { call file_size_gb as fq1_size { input: infile = fq1 } }
    #scatter (fq2 in fastq2.out) { call file_size_gb as fq2_size { input: infile = fq2 } }    
    call sum_fastq_size {input: R1 = fastq1.out, R2 = fastq2.out}
  
    # Split the fastq files into chunks for parallelization
    call chunk_fastq_files  { input:    sample_id = sample_id, 
                                        r1_in = fastq1.out, 
                                        r2_in = fastq2.out, 
                                        num_lines_per_chunk = 4 * num_reads_per_chunk, 
                                        disk_gb = 20 + sum_fastq_size.gb * 5 
                            }

    scatter (fastq_pair in chunk_fastq_files.fastq_pairs) {
          call microc_align {input: 
                        image_id = image_id, 
                        sample_id = sample_id, 
                        fastq_R1 = fastq_pair.left, 
                        fastq_R2 = fastq_pair.right,
                        sample_id = sample_id, 
                        reference_index = reference_bwa_idx, 
                        reference_index_prefix = reference_bwa_idx_prefix, 
                        chrom_sizes = chrom_sizes,
                        resource_monitor_script = resource_monitor_script,
                        top_monitor_script = top_monitor_script
                       }
      }

#     call merge_bams { input:         
#                              image_id = image_id, 
#                              sample_id = sample_id, 
#                              bams = microc_align.bam,
#                              resource_monitor_script = resource_monitor_script,
#                              disk_gb = 30 + sum_fastq_size.gb * 5 
#     }
    
    call merge_pairs {input: 
                             image_id = image_id, 
                             sample_id = sample_id, 
                             pairsams = microc_align.pairsam,
                             resource_monitor_script = resource_monitor_script,
                             disk_gb = 30 + sum_fastq_size.gb * 5 
    }
    
    call juicer_hic {input: 
                        image_id = image_id, 
                        sample_id = sample_id, 
                        chrom_sizes = chrom_sizes, 
                        mapped_pairs = merge_pairs.mapped_pairs
                    }

    call cooler {input: 
                        image_id = image_id, 
                        sample_id = sample_id, 
                        chrom_sizes = chrom_sizes, 
                        mapped_pairs = merge_pairs.mapped_pairs
                    }

    call version_info {input: image_id = image_id}

    call run_qc {input:
        image_id = image_id, 
        #mapped_pairs = merge_pairs.mapped_pairs,
        mapped_stats = merge_pairs.microc_stats,
        sample_id = sample_id
    }

    output {
        File stats = merge_pairs.microc_stats
        File mapped_pairs = merge_pairs.mapped_pairs
        File bam = merge_pairs.bam
        File bai = merge_pairs.bai
        File hic = juicer_hic.hic
        File raw_mcool = cooler.raw_mcool
        File balanced_mcool = cooler.balanced_mcool
        String pipeline_version = version_info.pipeline_version
        String reads_total = run_qc.reads_total
        String reads_mapped = run_qc.reads_mapped
        String reads_nodups = run_qc.reads_nodups
        String reads_cis_1kb = run_qc.reads_cis_1kb
        String reads_cis_10kb = run_qc.reads_cis_10kb
        Array[File] resources_align = microc_align.resources
        Array[File] top_align = microc_align.top
    }

}

task split_string_into_array {
    input {
        String? str
        String arr = "{ADDR[@]}"        
    }
    command {
        IFS=',' read -ra ADDR <<< "${str}"
        for i in "$${arr}"; do 
            echo "$i" | tr -d " "; 
        done        
    }
    runtime {
        docker: "ubuntu"
    }
    output {
        Array[String] out = read_lines(stdout())
    }
}

task chunk_fastq_files {
    input {
        String sample_id
        Array[File] r1_in
        Array[File] r2_in
        Int num_lines_per_chunk
        Int disk_gb
    }
    
    command {
        zcat -f ${sep=' ' r1_in} | split -d --suffix-length=3 -l ${num_lines_per_chunk} --additional-suffix='_R1.fastq' --filter='gzip > $FILE.gz' - ${sample_id}-
        zcat -f ${sep=' ' r2_in} | split -d --suffix-length=3 -l ${num_lines_per_chunk} --additional-suffix='_R2.fastq' --filter='gzip > $FILE.gz' - ${sample_id}-
    }
    
     runtime {
        continueOnReturnCode: false
        docker: "debian:stretch"
        cpu: 4
        disks: "local-disk " + disk_gb + " SSD"        
    }   
    output {
        Array[File] r1_out = glob("*_R1.fastq.gz")
        Array[File] r2_out = glob("*_R2.fastq.gz")        
        Array[Pair[File, File]] fastq_pairs = zip(r1_out, r2_out)
    }
}


# task file_size_gb {
#     input {File infile}
#     command {} 
#     runtime {
#         docker: "debian:stretch"
#         disks: "local-disk 250 SSD"
#     }
#     output {
#         Float gb = size(infile, "GB")
#     }      
# }


task sum_fastq_size {
    input {
        Array[File] R1
        Array[File] R2
        Int size_gb = round(size(R1, "GB") + size(R2, "GB"))
    }
    command <<<
        echo "Calculating fastq file size"
    >>>
    runtime {
        docker: "debian:stretch"
        disks: "local-disk 2000 SSD"
    }
    output {
        Int gb = size_gb
    }
}
          
task merge_fastqs {
    input {
        Array[File] fastq_r1
        Array[File] fastq_r2
        String memory = "20GB"
    }

    command {
        cat ~{sep=' ' fastq_r1} > R1.fastq.gz
        cat ~{sep=' ' fastq_r2} > R2.fastq.gz
    }

    runtime {
        docker: "ubuntu"
        cpu: 4
        memory: memory
        disks: "local-disk " + 300 + " SSD" 
    }

    output {
        File? fastq_out1 = "R1.fastq.gz"
        File? fastq_out2 = "R2.fastq.gz"
    }
}


task microc_align {
    input {
        String image_id
        String sample_id
        File fastq_R1
        File fastq_R2
        File? reference_index
        String? reference_index_prefix
        File chrom_sizes
        Int bwa_cores = 5
        String memory = "20GB"
        String disk = "100"
        String mapq = "20"
        String chunk_id = basename(fastq_R1, "_R1.fastq.gz")
        File resource_monitor_script
        File top_monitor_script
    }

    command {
        
        # Start process monitoring script (with top)
        chmod +x ${top_monitor_script}
        ${top_monitor_script} > ${chunk_id}.top.log &

        # Start resource monitoring script
        chmod +x ${resource_monitor_script}
        ${resource_monitor_script} > ${chunk_id}.resources.log &
        
        # Check if bwa index is provided as a tar.gz file ("reference_index")
        #  or a URI prefix ("reference_index_prefix", e.g. gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta.64)
        if [ ! -z "${reference_index}" ] && [ ! -z "${reference_index_prefix}" ]
        then
            echo "ERROR: Both reference_index and reference_index_prefix provided. Please provide only one."
            exit
        fi
        if [ ! -z "${reference_index}" ]
        then
            echo "Using provided reference .tar.gz: ${reference_index}"
            mkdir genome_index
            tar zxvf ${reference_index} -C genome_index
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
        BWT=$(find genome_index -name '*.bwt')
        GENOME_INDEX_FA="$(dirname $BWT)"/"$(basename $BWT .bwt)"
        echo "Using bwa index: $GENOME_INDEX_FA"
        
        bwa mem -5SP -T0 -t${bwa_cores} $GENOME_INDEX_FA ${fastq_R1} ${fastq_R2}| \
        pairtools parse --min-mapq ${mapq} --walks-policy 5unique \
        --max-inter-align-gap 30 --add-columns pos5,pos3,dist_to_5,dist_to_3,read_len \
        --nproc-in ${bwa_cores} --nproc-out ${bwa_cores} --chroms-path ${chrom_sizes} | \
        pairtools sort --nproc ${bwa_cores} > chunk.pairsam
    }

    runtime {
        docker: "us-central1-docker.pkg.dev/aryeelab/docker/microc:${image_id}"
        bootDiskSizeGb: 40
        cpu: bwa_cores
        memory: memory
        disks: "local-disk " + disk + " SSD"

    }

    output {
        File pairsam = "chunk.pairsam"
        File resources = "${chunk_id}.resources.log"
        File top = "${chunk_id}.top.log"
    }

}

task merge_bams {
    input {
        String image_id
        String sample_id
        Array[File] bams    
        File resource_monitor_script
        Int disk_gb
    }    
    command {
        chmod u+x ${resource_monitor_script}
        ${resource_monitor_script} > resources.log &

        samtools merge -o ${sample_id}.bam ${sep=' ' bams}
        samtools index ${sample_id}.bam
        
        
    }
    runtime {
        continueOnReturnCode: false
        docker: "us-central1-docker.pkg.dev/aryeelab/docker/microc:${image_id}"
        bootDiskSizeGb: 20
        cpu: 4            
        disks: "local-disk " + disk_gb + " SSD"        
    }
    output {
        File resources = "resources.log"
        File bam = "${sample_id}.bam"
        File bai = "${sample_id}.bam.bai"
    }
}

task merge_pairs {
    input {
        String image_id 
        String sample_id
        Array[File] pairsams    
        File resource_monitor_script
        Int disk_gb
    }
    command {
        chmod u+x ${resource_monitor_script}
        ${resource_monitor_script} > resources.log &

        pairtools merge --nproc 12 ${sep=' ' pairsams} | \
        pairtools dedup --nproc-in 6 --nproc-out 6 --mark-dups --output-stats ${sample_id}.stats.txt | \
        pairtools split --nproc-in 6 --nproc-out 6 --output-pairs ${sample_id}.mapped.pairs --output-sam -| \
        samtools view -bS -@6 | \
        samtools sort -@6 -o ${sample_id}.bam
                
        samtools index ${sample_id}.bam
        
    }
    runtime {
        continueOnReturnCode: false
        docker: "us-central1-docker.pkg.dev/aryeelab/docker/microc:${image_id}"
        bootDiskSizeGb: 20
        cpu: 16            
        disks: "local-disk " + disk_gb + " SSD"        
    }
    output {
        File resources = "resources.log"
        File mapped_pairs = "${sample_id}.mapped.pairs"
        File microc_stats = "${sample_id}.stats.txt"
        File bam = "${sample_id}.bam"
        File bai = "${sample_id}.bam.bai"       
    }
}


task juicer_hic {
    input {
        String image_id
        String sample_id
        File chrom_sizes
        File mapped_pairs
        Int cores = 2
        String memory = "40GB"
        String disk = "200"
    }

    command {
        java -Xmx32000m  -Djava.awt.headless=true -jar /usr/local/bin/juicer_tools_1.22.01.jar pre \
            --threads ${cores} \
            ${mapped_pairs} \
            ${sample_id}.hic \
            ${chrom_sizes}
    }

    runtime {
        docker: "us-central1-docker.pkg.dev/aryeelab/docker/juicer:${image_id}"
        bootDiskSizeGb: 40
        memory: memory
        disks: "local-disk " + disk + " SSD"
    }

    output {
        File hic = "${sample_id}.hic"
    
    }

}

task cooler {
    input {
        String image_id
        String sample_id
        File chrom_sizes
        File mapped_pairs
        Int resolution = "10000"
        String memory = "32"
        String disk = "2000"
    }

    command {
        cooler cload pairs -c1 2 -p1 3 -c2 4 -p2 5 ${chrom_sizes}:${resolution} ${mapped_pairs} ${sample_id}.cool
        cooler zoomify --resolutions ${resolution}N -o ${sample_id}.raw.mcool -p 4 ${sample_id}.cool
        cooler zoomify --resolutions ${resolution}N -o ${sample_id}.balanced.mcool -p 4 --balance --balance-args '--nproc 4' ${sample_id}.cool
    }

    runtime {
        docker: "us-central1-docker.pkg.dev/aryeelab/docker/cooler:${image_id}"
        cpu: 4
        memory: memory + "GB"
        disks: "local-disk " + disk + " SSD"
    }

    output {
        File raw_mcool = "${sample_id}.raw.mcool"   
        File balanced_mcool = "${sample_id}.balanced.mcool"         
    }
}

task version_info { 
    input {
        String image_id
    }
    
    command {
        cat /VERSION
    }
    
    runtime {
            continueOnReturnCode: false
            docker: "us-central1-docker.pkg.dev/aryeelab/docker/microc:${image_id}"
            cpu: 1
            memory: "1GB"
        }
    output {
        String pipeline_version = read_string(stdout())
    }
}

task run_qc {
    input {
        String image_id
        File mapped_stats
        String sample_id
    }

    command {
        cat ${mapped_stats} | grep -w "total" | cut -f2 
        cat ${mapped_stats} | grep -w "total_mapped" | cut -f2 
        cat ${mapped_stats} | grep -w "total_nodups" | cut -f2 
        cat ${mapped_stats} | grep -w "cis_1kb+" | cut -f2 
        cat ${mapped_stats} | grep -w "cis_10kb+" | cut -f2
    }

    runtime {
        docker: "us-central1-docker.pkg.dev/aryeelab/docker/utils:${image_id}"
        cpu: 1
        memory: "4GB"
        disks: "local-disk 10 SSD"
    }
    output {
        #File qc_stats_file = "${sample_id}_qc.zip"
        Array[String] qc_stats = read_lines(stdout())
        #String reads_total_with_commas = qc_stats[0]
        #String dist20kb_reads = qc_stats[1]
        #String dist20kb_percent = qc_stats[2]
        String reads_total = qc_stats[0]
        String reads_mapped = qc_stats[1]
        String reads_nodups = qc_stats[2]
        String reads_cis_1kb = qc_stats[3]
        String reads_cis_10kb = qc_stats[4]     
    }
}
