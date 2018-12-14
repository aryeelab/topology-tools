workflow merge_hic_samples {
    String set_id
    Array[File] hicpro_out_tars
    Array[Int] replicate_pairs
    String genome_id
    String genome_size
    String bin_size    
    
    File monitoring_script
    
    # Count total reads across samples/replicates
    call count_pairs { input: replicate_pairs = replicate_pairs }
     
    # Merge the HiC-Pro align results 
    call hicpro_merge { input: set_id = set_id, hicpro_out_tars = hicpro_out_tars, monitoring_script = monitoring_script, disk_gb = 5000}

    # Calculate the cis-long range percent metric
    call cis_long_range_percent {input: set_id = set_id, num_pairs = count_pairs.num_pairs, qc_stats = hicpro_merge.qc_stats}
    
    # Compute raw and ICE normalized hicpro contact matrices
    call hicpro_contact_matrices {input: set_id = set_id, all_valid_pairs = hicpro_merge.all_valid_pairs, genome_size = genome_size, bin_size=bin_size, monitoring_script = monitoring_script, disk_gb = 5000}

    # Generate Juicebox format .hic file
    call juicebox_hic {input: sample_id = set_id, all_valid_pairs = hicpro_merge.all_valid_pairs, genome_size = genome_size, monitoring_script = monitoring_script, disk_gb = 5000}

    # Generate sparseHiC format .rds file
    call sparseHic {input: sample_id = set_id, matrix_zip = hicpro_contact_matrices.matrix, genome_size = genome_size, genome_id = genome_id, monitoring_script = monitoring_script, disk_gb = 5000}
    
}

task count_pairs {
    Array[Int] replicate_pairs
        
    command <<<
        let "total_pairs=${sep='+' replicate_pairs}"
        echo $total_pairs
    >>>

    runtime {
        docker: "debian:stretch"
    }
    
    output {
        String num_pairs = read_string(stdout())
    }   
}


task hicpro_merge {
    String set_id
    Array[File] hicpro_out_tars
    
    File monitoring_script
    Int disk_gb
    
    command <<<
    
           chmod u+x ${monitoring_script}
            ${monitoring_script} > monitoring.log &
 
            mkdir -p hic_results/data/${set_id}
            for tar in ${sep=' ' hicpro_out_tars}; do
                tar xvf $tar --strip=3 --wildcards -C hic_results/data/${set_id} 'hic_results/data/*/*'
            done;
            
            # Run HiC-Pro
            /HiC-Pro/bin/HiC-Pro -s merge_persample -i hic_results/data -o . -c /HiC-Pro/config-hicpro.txt
    
            # Zip qc stats
            zip -j qc_stats.zip \
                bowtie_results/bwt2/${set_id}/${set_id}.mpairstat \
                hic_results/data/${set_id}/${set_id}_allValidPairs.mergestat \
                hic_results/data/${set_id}/${set_id}.mRSstat
            
            # Zip logs
            zip -j logs.zip logs/${set_id}/*

    >>>
    
            runtime {
            continueOnReturnCode: false
            docker: "aryeelab/hicpro:latest"
            cpu: 4            
            disks: "local-disk " + disk_gb + " HDD"        
        }
        
    output {
            File monitoring_log = "monitoring.log"
            File all_valid_pairs = "hic_results/data/${set_id}/${set_id}_allValidPairs"
            File qc_stats = "qc_stats.zip"
            File hicpro_logs = "logs.zip"
    }
}

task cis_long_range_percent {
    String set_id
    String num_pairs
    File qc_stats

    String dollar = "$"

        
    command <<<
        unzip -qq ${qc_stats}
        cis_long_range=${dollar}(cat ${set_id}_allValidPairs.mergestat | grep cis_longRange | cut -f2)
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
        String set_id
        
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

            mkdir -p pairs/${set_id}
            ln -sf   ${all_valid_pairs} pairs/${set_id}/
                        
            # Run HiC-Pro
            /HiC-Pro/bin/HiC-Pro -s build_contact_maps -s ice_norm -i pairs -o hicpro_out -c /HiC-Pro/config-hicpro.txt
            
            # Zip contact matrices
            zip -rj matrix.zip hicpro_out/hic_results/matrix/${set_id}/iced hicpro_out/hic_results/matrix/${set_id}/raw
                 
            # Zip logs
            zip -j logs.zip hicpro_out/logs/${set_id}/*

        >>>
                
        output {
            File matrix = "matrix.zip"
            File hicpro_logs = "logs.zip"
            File monitoring_log = "monitoring.log"
        }
        
        
        runtime {
            continueOnReturnCode: false
            docker: "aryeelab/hicpro:latest"
            memory: "32GB"
            disks: "local-disk " + disk_gb + " HDD"        
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
            disks: "local-disk " + disk_gb + " HDD"            
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
        Rscript /usr/local/bin/hicpro_to_sparsehic.R --sample_id=${sample_id}  --genome_size=/HiC-Pro/annotation/${genome_size} --genome_id=${genome_id} --cores=1

    >>>
    output {
        File rds = "${sample_id}.sparsehic.rds"
    }

    runtime {
            continueOnReturnCode: false
            docker: "aryeelab/hicpro:latest"
            memory: "104GB"
            disks: "local-disk " + disk_gb + " SSD"            
    }
}
