workflow preprocess_hic {
    String sample_id
    String r1_fastq
    String r2_fastq
    String genome_size
    File monitoring_script

    call split as fastq1 {input: str = r1_fastq}
    call split as fastq2 {input: str = r2_fastq}
    
    scatter (fq1 in fastq1.out) { call file_size_gb as fq1_size { input: infile = fq1 } }
    scatter (fq2 in fastq2.out) { call file_size_gb as fq2_size { input: infile = fq2 } }    
    call calculate_fastq_size {input: size1 = fq1_size.gb, size2 = fq2_size.gb}

    call count_pairs {input: r1_fastq = fastq1.out, disk_gb = 10 + calculate_fastq_size.gb * 3}
    call hicpro_align {input: sample_id = sample_id, r1_fastq = fastq1.out, r2_fastq = fastq2.out, genome_size = genome_size, monitoring_script = monitoring_script, disk_gb = 30 + calculate_fastq_size.gb * 10}
    call cis_long_range_percent {input: sample_id = sample_id, num_pairs = count_pairs.num_pairs, qc_stats = hicpro_align.qc_stats}
    call hicpro_contact_matrices {input: sample_id = sample_id, all_valid_pairs = hicpro_align.all_valid_pairs, genome_size = genome_size, monitoring_script = monitoring_script}
}

task split {
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

task file_size_gb {
  File infile  
  Float f = size(infile, "GB")
  String s = f 
  String string_before_decimal = sub(s, "\\..*", "") 
  Int final_int = string_before_decimal
  command {} 
  output {
        Int gb = final_int
  }      
}

task calculate_fastq_size {
    Array[Int] size1
    Array[Int] size2
    command { 
        let "gb=${sep=' + ' size1} + ${sep=' + ' size2}"
        echo $gb
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
        Array[File] r1_fastq
        Array[File] r2_fastq
        String sample_id
        
        File genome_index_tgz
        String genome_name
        String genome_fragment
        String genome_size
        String ligation_site
        
        String bowtie2_cores
        
        String read1_ext = "_R1_"
        String read2_ext = "_R2_" 
        String min_mapq="20"
        
        File monitoring_script
        
        String memory
        Int disk_gb
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

            # Set up input fastq directory 
            # Create symlinks with consistent naming to original fastqs
            mkdir -p fastq/${sample_id}
            r1_fastq=(${sep=' ' r1_fastq})
            r2_fastq=(${sep=' ' r2_fastq})            
            for i in ${dollar}{!r1_fastq[@]}; do
              # Read 1
              CMD="ln -s ${dollar}{r1_fastq[$i]} fastq/${sample_id}/${sample_id}_R1_$i.fastq.gz"
              echo "Symlinking fastq1: $CMD"
              $CMD
              # Read 2
              CMD="ln -s ${dollar}{r2_fastq[$i]} fastq/${sample_id}/${sample_id}_R2_$i.fastq.gz"
              echo "Symlinking fastq2: $CMD"
              $CMD
            done
                        
            # Run HiC-Pro
            /HiC-Pro/bin/HiC-Pro -s mapping -s proc_hic -s quality_checks -s merge_persample -i fastq -o hicpro_out -c /HiC-Pro/config-hicpro.txt
            
            # Zip qc stats
            zip -j qc_stats.zip \
                hicpro_out/hic_results/data/${sample_id}/${sample_id}_allValidPairs.mergestat \
                hicpro_out/hic_results/data/${sample_id}/${sample_id}.mRSstat
            
            # Zip logs
            zip -j logs.zip hicpro_out/logs/${sample_id}/*

        >>>
                
        output {
            File all_valid_pairs = "hicpro_out/hic_results/data/${sample_id}/${sample_id}_allValidPairs"
            File qc_stats = "qc_stats.zip"
            File hicpro_logs = "logs.zip"
            File monitoring_log = "monitoring.log"
        }
                
        runtime {
            continueOnReturnCode: false
            docker: "aryeelab/hicpro:latest"
            cpu: cpu
            memory: memory
            disks: "local-disk " + disk_gb + " SSD"        
            preemptible: preemptible
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
            memory: "8GB"
            disks: "local-disk 200 SSD"
        }
}




