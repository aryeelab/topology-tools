version 1.0


workflow microc {
	String pipeline_ver = 'v0.0.0'

	meta {
		version: 'v0.0.0'

		author: ''
		email: ''
		description: 'Process microC data. Inspired in pipeline from: https://micro-c.readthedocs.io/en/latest/index.html'
		organization: ''

		default_docker: 'salvacasani/microc:latest'

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
		String fastq_R1
		String fastq_R2
		String reference_bwa_idx
		String chroms_path
		String docker = 'salvacasani/microc:latest'
	}

	call split_string_into_array as fastq1 {input : str = fastq_R1}
	call split_string_into_array as fastq2 {input : str = fastq_R2}

	call merge_fastqs {input : fastq_r1 = fastq1.out, fastq_r2 = fastq2.out}

	call microc_align {input : sample_id = sample_id, fastq_R1 = merge_fastqs.fastq_out1, fastq_R2 = merge_fastqs.fastq_out2,
		sample_id = sample_id, reference_index = reference_bwa_idx, chroms_path = chroms_path, docker = docker
	}

	output {
		File stats = microc_align.microc_stats
		File mapped_pairs = microc_align.mapped_pairs
		File bam = microc_align.bam
	}

}

task split_string_into_array {
    input {
        String str
        String arr = "{ADDR[@]}"        
    }
    command {
        IFS=';' read -ra ADDR <<< "${str}"
        for i in "$${arr}"; do echo "$i"; done        
    }
    runtime {
        docker: "ubuntu"
    }
    output {
        Array[String] out = read_lines(stdout())
    }
}

task merge_fastqs {
	input {
		Array[File] fastq_r1
		Array[File] fastq_r2
	}

	command {
		zcat -f ~{sep=' ' fastq_r1} | gzip > R1.fastq.gz
		zcat -f ~{sep=' ' fastq_r2} | gzip > R2.fastq.gz
	}

	runtime {
		docker: "salvacasani/microc:latest"
		cpu: 4
		disks: "local-disk " + 30 + " SSD" 
	}

	output {
		File fastq_out1 = "R1.fastq.gz"
		File fastq_out2 = "R2.fastq.gz"
	}
}


task microc_align {
	input {
		String sample_id
		File fastq_R1
		File fastq_R2
		File reference_index
		File chroms_path
		Int bwa_cores = 5
		String genome = "hg38"
		String docker
	}

	command {
		mkdir /cromwell_root/genome_index
		tar zxvf ${reference_index} -C /cromwell_root/genome_index
		bwa mem -5SP -T0 -t ${bwa_cores} /cromwell_root/genome_index/bwa_indexes/${genome} ${fastq_R1} ${fastq_R2}| \
		pairtools parse --min-mapq 40 --walks-policy 5unique \
		--max-inter-align-gap 30 --nproc-in ${bwa_cores} --nproc-out ${bwa_cores} --chroms-path ${chroms_path} | \
		pairtools sort --nproc ${bwa_cores}|pairtools dedup --nproc-in ${bwa_cores} \
		--nproc-out ${bwa_cores} --mark-dups --output-stats stats.txt|pairtools split --nproc-in ${bwa_cores} \
		--nproc-out ${bwa_cores} --output-pairs mapped.pairs --output-sam -|samtools view -bS -@${bwa_cores} | \
		samtools sort -@${bwa_cores} -o ${sample_id}.bam
	}

	runtime {
		docker: docker
		bootDiskSizeGb: 20
		cpu: bwa_cores
		disks: "local-disk 40 SSD"

	}

	output {
		File microc_stats = read_lines("stats.txt")
		File mapped_pairs = read_lines("mapped.pairs")
		File bam = "${sample_id}.bam"
	}

}


