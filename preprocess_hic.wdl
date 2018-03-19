workflow preprocess_hic {
    String sample_id
    String r1_fastq
    String r2_fastq
    Int num_reads_per_chunk
    String genome_id
    String genome_size
    String bin_size
    
    File monitoring_script

    # Split the comma-separated string of fastq file names into an array
    call split_string_into_array as fastq1 {input: str = r1_fastq}
    call split_string_into_array as fastq2 {input: str = r2_fastq}
        
    # Calculate the total fastq file size
    scatter (fq1 in fastq1.out) { call file_size_gb as fq1_size { input: infile = fq1 } }
    scatter (fq2 in fastq2.out) { call file_size_gb as fq2_size { input: infile = fq2 } }    
    call sum_fastq_size {input: size1 = fq1_size.gb, size2 = fq2_size.gb}

    # Count the number of read pairs
    call count_pairs {input: r1_fastq = fastq1.out, disk_gb = 10 + sum_fastq_size.gb * 3}

    # Split the fastq files into chunks for parallelization
    call split_fastq_files  { input: sample_id = sample_id, r1_in = fastq1.out, r2_in = fastq2.out, num_lines_per_chunk = 4 * num_reads_per_chunk, num_pairs = count_pairs.num_pairs, disk_gb = 20 + sum_fastq_size.gb * 3 }

    # Run HiC-Pro align on each fastq chunk
    scatter (fastq_pair in split_fastq_files.fastq_pairs) { call hicpro_align {input: sample_id = sample_id, r1_fastq = fastq_pair.left, r2_fastq = fastq_pair.right, genome_size = genome_size, monitoring_script = monitoring_script} }
    
    # Merge the HiC-Pro align results 
    call hicpro_merge { input: sample_id = sample_id, hicpro_out_tars = hicpro_align.hicpro_out, monitoring_script = monitoring_script, disk_gb = 30 + sum_fastq_size.gb * 10}

    # Calculate the cis-long range percent metric
    call cis_long_range_percent {input: sample_id = sample_id, num_pairs = count_pairs.num_pairs, qc_stats = hicpro_merge.qc_stats}
    
    # Compute raw and ICE normalized hicpro contact matrices
    call hicpro_contact_matrices {input: sample_id = sample_id, all_valid_pairs = hicpro_merge.all_valid_pairs, genome_size = genome_size, bin_size=bin_size, monitoring_script = monitoring_script, disk_gb = 30 + sum_fastq_size.gb * 10}

    # Generate Juicebox format .hic file
    call juicebox_hic {input: sample_id = sample_id, all_valid_pairs = hicpro_merge.all_valid_pairs, genome_size = genome_size, monitoring_script = monitoring_script, disk_gb = 30 + sum_fastq_size.gb}

    # Generate sparseHiC format .rds file
    call sparseHic {input: sample_id = sample_id, matrix_zip = hicpro_contact_matrices.matrix, genome_size = genome_size, genome_id = genome_id, monitoring_script = monitoring_script, disk_gb = 30 + sum_fastq_size.gb}

    # Generate balanced and unbalanced cooler files
    call cooler {input: sample_id = sample_id, all_valid_pairs = hicpro_merge.all_valid_pairs, genome_size = genome_size, bin_size=bin_size, monitoring_script = monitoring_script, disk_gb = 30 + sum_fastq_size.gb * 3}
}

task split_string_into_array {
    String str 
    String arr = "{ARR[@]}"
    command <<<
        IFS=',' read -ra ARR <<< "${str}"
        for i in "$${arr}"; do
            echo "$i" | tr -d " " >> out
        done
    >>>

    runtime {
        docker: "debian:stretch"
    }
    
    output {
        Array[String] out = read_lines("out")
    }
}

task split_fastq_files {
    String sample_id
    Array[File] r1_in
    Array[File] r2_in
    Int num_lines_per_chunk
    Int num_pairs
    Int disk_gb
    
    command {
        zcat ${sep=' ' r1_in} | split -d --suffix-length=3 -l ${num_lines_per_chunk} --additional-suffix='_R1.fastq' --filter='gzip > $FILE.gz' - ${sample_id}-
        zcat ${sep=' ' r2_in} | split -d --suffix-length=3 -l ${num_lines_per_chunk} --additional-suffix='_R2.fastq' --filter='gzip > $FILE.gz' - ${sample_id}-
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

task file_size_gb {
    File infile  
    command {} 
    runtime {
        docker: "debian:stretch"
        disks: "local-disk 100 SSD"
    }
    output {
        Float gb = size(infile, "GB")
    }      
}

task sum_fastq_size {
    Array[Float] size1
    Array[Float] size2
    command { 
        echo "(${sep=' + ' size1} + ${sep=' + ' size2})/1" | sed -e 's/[eE]+*/\*10\^/g' | bc
    }
    runtime {
        docker: "aryeelab/hicpro:latest"
    }
    output {
        Int gb = read_int(stdout())        
    }      
}

task count_pairs {
    Array[File] r1_fastq
    Int disk_gb
    String dollar = "$"
        
    command <<<
        num_lines=`zcat ${sep=' ' r1_fastq} | wc -l`
        let "num_pairs=${dollar}num_lines/4"
        echo ${dollar}num_pairs
    >>>

    runtime {
        docker: "debian:stretch"
        disks: "local-disk " + disk_gb + " SSD"
    }
    
    output {
        Int num_pairs = read_int(stdout())
    }   
}

task hicpro_align {
        File r1_fastq
        File r2_fastq
        String sample_id
        
        File genome_index_tgz
        String genome_name
        String genome_fragment
        String genome_size
        String ligation_site
        
        String bowtie2_cores
        
        String read1_ext = "_R1"
        String read2_ext = "_R2" 
        String min_mapq="20"
        
        File monitoring_script
        
        String memory
        Int cpu
        Int preemptible
        
        String dollar = "$"
        String at = "@"
        
        command <<<
                
            chmod u+x ${monitoring_script}
            ${monitoring_script} > monitoring.log &

            mkdir $PWD/bowtie2_index
            tar zxvf ${genome_index_tgz} -C $PWD/bowtie2_index

            # Set up hicpro config file
            CONFIG=/HiC-Pro/config-hicpro.txt
            sed -i "s|BOWTIE2_IDX_PATH.*|BOWTIE2_IDX_PATH = $PWD/bowtie2_index|" $CONFIG
            sed -i "s/N_CPU.*/N_CPU = ${bowtie2_cores}/" $CONFIG
            sed -i "s/PAIR1_EXT.*/PAIR1_EXT = ${read1_ext}/" $CONFIG
            sed -i "s/PAIR2_EXT.*/PAIR2_EXT = ${read2_ext}/" $CONFIG
            sed -i "s/MIN_MAPQ.*/MIN_MAPQ = ${min_mapq}/" $CONFIG
            sed -i "s/GENOME_FRAGMENT.*/GENOME_FRAGMENT = ${genome_fragment}/" $CONFIG
            sed -i "s/LIGATION_SITE.*/LIGATION_SITE = ${ligation_site}/" $CONFIG
            sed -i "s/REFERENCE_GENOME.*/REFERENCE_GENOME = ${genome_name}/" $CONFIG
            sed -i "s/GENOME_SIZE.*/GENOME_SIZE = ${genome_size}/" $CONFIG

            # Set up input fastq directory using symlinks to fastqs
            mkdir -p rawdata/${sample_id}
            ln -s ${r1_fastq} rawdata/${sample_id}/
            ln -s ${r2_fastq} rawdata/${sample_id}/
            
            mkdir tmp
            mkdir logs     
            #echo ${sample_id}/`basename ${r1_fastq}` > inputfiles.txt
            (cd rawdata && ls */*_R1.fastq.gz) > inputfiles.txt
            export FASTQFILE=inputfiles.txt
            export LSB_JOBINDEX=1
    
            # Run HiC-Pro
            # all_sub : configure bowtie_global bowtie_local bowtie_combine mapping_stat bowtie_pairing mapped_2hic_fragments 
            make --file //HiC-Pro/scripts/Makefile CONFIG_FILE=/HiC-Pro/config-hicpro.txt CONFIG_SYS=//HiC-Pro/config-system.txt all_sub 2>&1
            
            #partname=$(echo $(basename imr90-rep1/part-01_R1.fastq.gz) | sed 's/_R1.fastq.gz//')
            tar -cvpf hicpro_out.tar bowtie_results/bwt2 hic_results logs
        >>>
                
        output {
            File hicpro_out = "hicpro_out.tar"
            File monitoring_log = "monitoring.log"
        }
                
        runtime {
            continueOnReturnCode: false
            docker: "aryeelab/hicpro:latest"
            cpu: cpu
            memory: memory
            disks: "local-disk 20 SSD"        
            preemptible: preemptible
        }
}


task hicpro_merge {
    String sample_id
    Array[File] hicpro_out_tars
    
    File monitoring_script
    Int disk_gb
    
    command <<<
    
           chmod u+x ${monitoring_script}
            ${monitoring_script} > monitoring.log &
 
            for tar in ${sep=' ' hicpro_out_tars}; do
                tar xvf $tar
            done;
            
            # Run HiC-Pro
            /HiC-Pro/bin/HiC-Pro -s merge_persample -i hic_results/data -o . -c /HiC-Pro/config-hicpro.txt
    
            # Zip qc stats
            zip -j qc_stats.zip \
                bowtie_results/bwt2/${sample_id}/${sample_id}.mpairstat \
                hic_results/data/${sample_id}/${sample_id}_allValidPairs.mergestat \
                hic_results/data/${sample_id}/${sample_id}.mRSstat
            
            # Zip logs
            zip -j logs.zip logs/${sample_id}/*

    >>>
    
            runtime {
            continueOnReturnCode: false
            docker: "aryeelab/hicpro:latest"
            cpu: 4            
            disks: "local-disk " + disk_gb + " SSD"        
        }
        
    output {
            File monitoring_log = "monitoring.log"
            File all_valid_pairs = "hic_results/data/${sample_id}/${sample_id}_allValidPairs"
            File qc_stats = "qc_stats.zip"
            File hicpro_logs = "logs.zip"
    }
}

task cis_long_range_percent {
    String sample_id
    Int num_pairs
    File qc_stats

    String dollar = "$"

        
    command <<<
        unzip -qq ${qc_stats}
        cis_long_range=${dollar}(cat ${sample_id}_allValidPairs.mergestat | grep cis_longRange | cut -f2)
        let "percent=100*$cis_long_range/${num_pairs}"
        echo $percent
    >>>

    runtime {
        docker: "aryeelab/hicpro:latest"
    }
    
    output {
        Int percent = read_int(stdout())
    }  
}

task hicpro_contact_matrices {
        File all_valid_pairs
        String sample_id
        
        String genome_size
        String bin_size

        Int disk_gb
        File monitoring_script
                        
        command <<<
                
            chmod u+x ${monitoring_script}
            ${monitoring_script} > monitoring.log &

            CONFIG=/HiC-Pro/config-hicpro.txt
            sed -i "s/BIN_SIZE.*/BIN_SIZE = ${bin_size}/" $CONFIG
            sed -i "s/GENOME_SIZE.*/GENOME_SIZE = ${genome_size}/" $CONFIG

            mkdir -p pairs/${sample_id}
            ln -sf   ${all_valid_pairs} pairs/${sample_id}/
                        
            # Run HiC-Pro
            /HiC-Pro/bin/HiC-Pro -s build_contact_maps -s ice_norm -i pairs -o hicpro_out -c /HiC-Pro/config-hicpro.txt
            
            # Zip contact matrices
            zip -rj matrix.zip hicpro_out/hic_results/matrix/${sample_id}/iced hicpro_out/hic_results/matrix/${sample_id}/raw
                 
            # Zip logs
            zip -j logs.zip hicpro_out/logs/${sample_id}/*

        >>>
                
        output {
            File matrix = "matrix.zip"
            File hicpro_logs = "logs.zip"
            File monitoring_log = "monitoring.log"
        }
        
        
        runtime {
            continueOnReturnCode: false
            docker: "aryeelab/hicpro:latest"
            memory: "208GB"
            disks: "local-disk " + disk_gb + " SSD"        
        }
}

task juicebox_hic {
    String sample_id
    File all_valid_pairs
    String genome_size
    
    Int disk_gb
    File monitoring_script
    
    command <<<
        chmod u+x ${monitoring_script}
        ${monitoring_script} > monitoring.log &

        /HiC-Pro/bin/utils/hicpro2juicebox.sh -i ${all_valid_pairs} -g /HiC-Pro/annotation/${genome_size} -j /usr/local/juicer/juicer_tools.1.7.6_jcuda.0.8.jar
        # Rename output .hic file
        mv ${sample_id}_allValidPairs.hic ${sample_id}.hic
    >>>
    output {
        File juicebox_hic = "${sample_id}.hic"
    }

    runtime {
            continueOnReturnCode: false
            docker: "aryeelab/hicpro:latest"
            memory: "60GB"
            disks: "local-disk " + disk_gb + " SSD"            
    }
}

task sparseHic {
    String sample_id
    String genome_size
    String genome_id
    File matrix_zip
    
    Int disk_gb
    File monitoring_script
    
    command <<<
        chmod u+x ${monitoring_script}           
        ${monitoring_script} > monitoring.log &
    
        unzip ${matrix_zip}
        Rscript /usr/local/bin/hicpro_to_sparsehic.R --sample_id=${sample_id}  --genome_size=/HiC-Pro/annotation/${genome_size} --genome_id=${genome_id} --cores=2

    >>>
    output {
        File rds = "${sample_id}.sparsehic.rds"
    }

    runtime {
            continueOnReturnCode: false
            docker: "aryeelab/hicpro:latest"
            memory: "60GB"
            disks: "local-disk " + disk_gb + " SSD"            
    }
}

task cooler {
    String sample_id
    File all_valid_pairs
    String genome_size
    String bin_size
    
    Int disk_gb  
    File monitoring_script

    command <<<

        chmod u+x ${monitoring_script}
        ${monitoring_script} > monitoring.log &
    
        echo`date`: Choosing smallest bin size of ${bin_size}
        RES=$(echo "${bin_size}" | tr " " "\n" | sort | head -n1)
    
        echo `date`: Starting makebins
        cooler makebins /annotation/${genome_size} $RES > bins.bed

        echo `date`: Starting cooler csort
        cooler csort --nproc 3 --chrom1 2 --pos1 3 --chrom2 5 --pos2 6 \
             -o allValidPairs.sorted  ${all_valid_pairs} /annotation/${genome_size}

        echo `date`: Starting cooler cload pairix
        cooler cload pairix bins.bed allValidPairs.sorted ${sample_id}.cool

        echo "`date`: Starting cooler zoomify (unbalanced)"
        cooler zoomify --no-balance --out ${sample_id}.mcool ${sample_id}.cool

        echo "`date`: Starting cooler zoomify (balanced)"
        cooler zoomify --balance --out ${sample_id}.balanced.mcool ${sample_id}.cool

        echo `date`: Done
    >>>

    runtime {
        continueOnReturnCode: false    
        docker: "aryeelab/cooler:latest"
        memory: "60GB"
        disks: "local-disk " + disk_gb + " SSD"        

    }
    
    output {
        File mcool = "${sample_id}.mcool"
        File monitoring_log = "monitoring.log"
    }
}
