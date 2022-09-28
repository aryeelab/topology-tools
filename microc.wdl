version 1.0


workflow microc {
	String pipeline_ver = 'v0.0.0'

	meta {
		version: 'v0.0.0'

		author: ''
		email: ''
		description: 'Process microC data. Inspired in pipeline from: https://micro-c.readthedocs.io/en/latest/index.html'
		organization: ''

		default_docker: 'us-central1-docker.pkg.dev/aryeelab/docker/microc:latest'

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
		String docker = 'us-central1-docker.pkg.dev/aryeelab/docker/microc:latest'
	}

	call split_string_into_array as fastq1 {input : str = fastq_R1}
	call split_string_into_array as fastq2 {input : str = fastq_R2}

	call merge_fastqs {input : fastq_r1 = fastq1.out, fastq_r2 = fastq2.out}

	call microc_align {input : sample_id = sample_id, fastq_R1 = merge_fastqs.fastq_out1, fastq_R2 = merge_fastqs.fastq_out2,
		sample_id = sample_id, reference_index = reference_bwa_idx, chroms_path = chroms_path, docker = docker
	}
	
	call juicer_hic {input : sample_id = sample_id, chroms_path = chroms_path, mapped_pairs = microc_align.mapped_pairs}

	output {
		File stats = microc_align.microc_stats
		File mapped_pairs = microc_align.mapped_pairs
		File bam = microc_align.bam
		File bai = microc_align.bai
		File hic = juicer_hic.hic
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
		String memory = "20GB"
	}

	command {
		zcat -f ~{sep=' ' fastq_r1} | gzip > R1.fastq.gz
		zcat -f ~{sep=' ' fastq_r2} | gzip > R2.fastq.gz
	}

	runtime {
		docker: "ubuntu"
		cpu: 4
		memory: memory
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
		String docker
		String memory = "20GB"
	}

	command {
		mkdir genome_index
		tar zxvf ${reference_index} -C genome_index
		# Get genome index fasta name. 
		BWT=$(find genome_index -name '*.bwt')	
		GENOME_INDEX_FA="$(dirname $BWT)"/"$(basename $BWT .bwt)"
		echo "Using bwa index: $GENOME_INDEX_FA"
		
		
		bwa mem -5SP -T0 -t${bwa_cores} $GENOME_INDEX_FA ${fastq_R1} ${fastq_R2}| \
		pairtools parse --min-mapq 40 --walks-policy 5unique \
		--max-inter-align-gap 30 --nproc-in ${bwa_cores} --nproc-out ${bwa_cores} --chroms-path ${chroms_path} | \
		pairtools sort --nproc ${bwa_cores} | pairtools dedup --nproc-in ${bwa_cores} \
		--nproc-out ${bwa_cores} --mark-dups --output-stats stats.txt | pairtools split --nproc-in ${bwa_cores} \
		--nproc-out ${bwa_cores} --output-pairs mapped.pairs --output-sam -|samtools view -bS -@${bwa_cores} | \
		samtools sort -@${bwa_cores} -o ${sample_id}.bam; samtools index ${sample_id}.bam
	}

	runtime {
		docker: docker
		bootDiskSizeGb: 40
		cpu: bwa_cores
		memory: memory
		disks: "local-disk 60 SSD"

	}

	output {
		File microc_stats = "stats.txt"
		File mapped_pairs = "mapped.pairs"
		File bam = "${sample_id}.bam"
		File bai = "${sample_id}.bam.bai"
	}

}

task juicer_hic {
	input {
		String sample_id
		File chroms_path
		File mapped_pairs
		Int cores = 2
		String memory = "40GB"
	}

	command {
		java -Xmx32000m  -Djava.awt.headless=true -jar /usr/local/bin/juicer_tools_1.22.01.jar pre \
			--threads ${cores} \
			${mapped_pairs} \
			${sample_id}.hic \
			${chroms_path}
	}

	runtime {
		docker: "us-central1-docker.pkg.dev/aryeelab/docker/juicer:latest"
		bootDiskSizeGb: 40
		memory: memory
		disks: "local-disk 200 SSD"
	}

	output {
		File hic = "${sample_id}.hic"
	
	}

}
