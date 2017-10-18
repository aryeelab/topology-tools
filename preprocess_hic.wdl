workflow preprocess_hic {
        call hicpro { }
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
        
        String read1_ext
        String read2_ext
        String bin_size="1000 2000"
        String min_mapq="20"
        
        File monitoring_script
        
        String memory
        String disks
        Int cpu
        Int preemptible
        
        command {
                
            chmod u+x ${monitoring_script}
            #${monitoring_script} > monitoring.log &

            mkdir /bowtie2_index
            tar zxvf ${genome_index_tgz} -C /bowtie2_index

            # Set up hicpro config file
            CONFIG=/HiC-Pro/config-hicpro.txt
            sed -i "s/BOWTIE2_IDX_PATH.*/BOWTIE2_IDX_PATH = \/bowtie2_index/" $CONFIG
            sed -i "s/N_CPU.*/N_CPU = ${cpu}/" $CONFIG
            sed -i "s/PAIR1_EXT.*/PAIR1_EXT = ${read1_ext}/" $CONFIG
            sed -i "s/PAIR2_EXT.*/PAIR2_EXT = ${read2_ext}/" $CONFIG
            sed -i "s/BIN_SIZE.*/BIN_SIZE = ${bin_size}/" $CONFIG
            sed -i "s/MIN_MAPQ.*/MIN_MAPQ = ${min_mapq}/" $CONFIG
            sed -i "s/GENOME_FRAGMENT.*/GENOME_FRAGMENT = ${genome_fragment}/" $CONFIG
            sed -i "s/LIGATION_SITE.*/LIGATION_SITE = ${ligation_site}/" $CONFIG
            sed -i "s/REFERENCE_GENOME.*/REFERENCE_GENOME = ${genome_name}/" $CONFIG
            sed -i "s/GENOME_SIZE.*/GENOME_SIZE = ${genome_size}/" $CONFIG

            # Set up input fastq directory
            mkdir -p /fastq/${sample_id}
            for fq in ${sep=' ' r1_fastq}; do ln -s $fq /fastq/${sample_id}/; done
            for fq in ${sep=' ' r2_fastq}; do ln -s $fq /fastq/${sample_id}/; done

            # Run HiC-Pro
            /HiC-Pro/bin/HiC-Pro -i /fastq -o hicpro_out -c /HiC-Pro/config-hicpro.txt
            
            # Zip contact matrices
            zip -rj matrix.zip hicpro_out/hic_results/matrix/${sample_id}/iced hicpro_out/hic_results/matrix/${sample_id}/raw
            
            # Zip logs
            cd hicpro_out/logs/${sample_id}
            zip logs.zip *

        }
                
        output {
            File matrix = "matrix.zip"
            File logs = "hicpro_out/logs/${sample_id}/logs.zip"
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
