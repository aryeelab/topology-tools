workflow preprocess_hic {
    String sample_id
    String r1_fastq
    String r2_fastq
    call split as fastq1 {input: str = r1_fastq}
    call split as fastq2 {input: str = r2_fastq}
    call hicpro {input: sample_id = sample_id, r1_fastq = fastq1.out, r2_fastq = fastq2.out}
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

task hicpro {
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
        String bin_size
        String min_mapq="20"
        
        File monitoring_script
        
        String memory
        String disks
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
            sed -i "s/BIN_SIZE.*/BIN_SIZE = ${bin_size}/" $CONFIG
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
            /HiC-Pro/bin/HiC-Pro -i fastq -o hicpro_out -c /HiC-Pro/config-hicpro.txt
            
            # Zip contact matrices
            zip -rj matrix.zip hicpro_out/hic_results/matrix/${sample_id}/iced hicpro_out/hic_results/matrix/${sample_id}/raw
            
            # Zip qc plots
            zip -rj qc_plots.zip hicpro_out/hic_results/pic/${sample_id}
            
            # Zip logs
            cd hicpro_out/logs/${sample_id}
            zip logs.zip *

        >>>
                
        output {
            File matrix = "matrix.zip"
            File all_valid_pairs = "hicpro_out/hic_results/data/${sample_id}/${sample_id}_allValidPairs"
            File hicpro_logs = "hicpro_out/logs/${sample_id}/logs.zip"
            File monitoring_log = "monitoring.log"
        }
        
        runtime {
            continueOnReturnCode: false
            docker: "aryeelab/hicpro:latest"
            memory: memory
            disks: disks
            cpu: cpu
            preemptible: preemptible
        }
}
