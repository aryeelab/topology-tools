workflow microc {
	String pipeline_ver = 'v0.0.0'

	meta {
        version: 'v0.0.0'

        author: ''
        email: ''
        description: 'Process microC data. Inspired in pipeline from: https://micro-c.readthedocs.io/en/latest/index.html'
        organization: ''

        specification_document: ''

        default_docker: 'salvacasani/microc:latest'
        default_singularity: ''
        default_conda: ''
        croo_out_def: ''

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
                help: 'The alignment is done using bwa, and requires the reference genome in fasta format, indexed. The number of cores is hardoded for now.'
            }
        }
    }

    input {
		String sample_id

    	String fastq_R1
    	String fastq_R2

        String ref_genome='hg38.fa'

        String docker = 'salvacasani/microc:latest'
    }

    RuntimeEnvironment runtime_environment = {
        'docker': docker
    }

    call split_string_into_array as fastq1 {input: str = fastq_R1}
    call split_string_into_array as fastq2 {input: str = fastq_R2}

    call merge_fastqs {input: fastq_R1 = fastq1.out, fastq_R2 = fastq2.out}

    call microc {input: sample_id = sample_id, fastq_R1 = merge_fastqs.fastq_out1, fastq_R2 = merge_fastqs.fastq_out2}

    output {
    	File microc.microc_stats
    	File microc.mapped_pairs
    	File microc.bam
    }

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

task merge_fastqs {
	Array[File] fastq_r1
	Array[File] fastq_r2
	String sample_id

	command {
		zcat -f ${sep=' ' fastq_r1} | gzip > '${sample_id}_R1.fastq.gz'
		zcat -f ${sep=' ' fastq_r2} | gzip > '${sample_id}_R2.fastq.gz'
	}

	runtime {
        docker: "debian:stretch"
        cpu: 4
        disks: "local-disk " + 30 + " SSD" 
    }
    
    output {
        File fastq_out1 = '${sample_id}_R1.fastq.gz'
        File fastq_out2 = '${sample_id}_R2.fastq.gz'
    }
}


task microc {
	String sample_id
	File fastq_R1
	File fastq_R2
	File reference_fa
	File reference_genome
	Int bwa_cores = 5

    command {
        bwa mem -5SP -T0 -t bwa_cores reference_fa fastq_R1 fastq_R2| \
        pairtools parse --min-mapq 40 --walks-policy 5unique \
        --max-inter-align-gap 30 --nproc-in bwa_cores --nproc-out bwa_cores --chroms-path reference_genome | \
        pairtools sort --nproc bwa_cores|pairtools dedup --nproc-in bwa_cores \
        --nproc-out bwa_cores --mark-dups --output-stats stats.txt|pairtools split --nproc-in bwa_cores \
        --nproc-out bwa_cores --output-pairs mapped.pairs --output-sam -|samtools view -bS -@bwa_cores | \
        samtools sort -@mapped.pairs -o '${sample_id}.bam';samtools index mapped.PT.bam
    }

    runtime {
        docker: runtime_environment.docker
        bootDiskSizeGb: 20
        cpu: bwa_cores
        disks: "local-disk 40 SSD"
        preemptible: preemptible

    }
    
    output {
    	File microc_stats = stats.txt
    	File mapped_pairs = mapped.pairs
    	File bam = '${sample_id}.bam'
    }

}


