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
		File? fastq_R1
		File? fastq_R2
		File? reference_bwa_idx
		String? reference_bwa_idx_prefix
		File chrom_sizes
		Boolean merge = true
	}

	if ( merge ) {
		call split_string_into_array as fastq1 {input: 
												str = fastq_R1
											}
		call split_string_into_array as fastq2 {input: 
												str = fastq_R2
											}
		call merge_fastqs {input: 
							fastq_r1 = fastq1.out, 
							fastq_r2 = fastq2.out}
	}
	
	File? fastq_R1_align = if merge then merge_fastqs.fastq_out1
					else fastq_R1
	File? fastq_R2_align = if merge then merge_fastqs.fastq_out2
					else fastq_R2

	call microc_align {input: 
						image_id = image_id, 
						sample_id = sample_id, 
						fastq_R1 = fastq_R1_align, 
						fastq_R2 = fastq_R2_align,
						sample_id = sample_id, 
						reference_index = reference_bwa_idx, 
						reference_index_prefix = reference_bwa_idx_prefix, 
						chrom_sizes = chrom_sizes
					   }
	
	call juicer_hic {input: 
						image_id = image_id, 
						sample_id = sample_id, 
						chrom_sizes = chrom_sizes, 
						mapped_pairs = microc_align.mapped_pairs
					}

	call version_info {input: image_id = image_id}

	call run_qc {input:
		image_id = image_id, 
		mapped_pairs = microc_align.mapped_pairs,
		mapped_stats = microc_align.microc_stats,
		sample_id = sample_id
	}

	output {
		File stats = microc_align.microc_stats
		File mapped_pairs = microc_align.mapped_pairs
		File bam = microc_align.bam
		File bai = microc_align.bai
		File hic = juicer_hic.hic
		String pipeline_version = version_info.pipeline_version
		File qcstats = run_qc.qc_stats_file
		String perc_20kb = run_qc.dist20kb
	}

}

task split_string_into_array {
    input {
        String? str
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
		File? fastq_out1 = "R1.fastq.gz"
		File? fastq_out2 = "R2.fastq.gz"
	}
}


task microc_align {
	input {
		String image_id
		String sample_id
		File? fastq_R1
		File? fastq_R2
		File? reference_index
		String? reference_index_prefix
		File chrom_sizes
		Int bwa_cores = 5
		String memory = "20GB"
		String disk = "500"
		String mapq = "20"
	}

	command {
		
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
		--max-inter-align-gap 30 --nproc-in ${bwa_cores} --nproc-out ${bwa_cores} --chroms-path ${chrom_sizes} | \
		pairtools sort --nproc ${bwa_cores} | pairtools dedup --backend cython --nproc-in ${bwa_cores} \
		--nproc-out ${bwa_cores} --mark-dups --output-stats stats.txt | pairtools split --nproc-in ${bwa_cores} \
		--nproc-out ${bwa_cores} --output-pairs mapped.pairs --output-sam -|samtools view -bS -@${bwa_cores} | \
		samtools sort -@${bwa_cores} -o ${sample_id}.bam; samtools index ${sample_id}.bam
	}

	runtime {
		docker: "us-central1-docker.pkg.dev/aryeelab/docker/microc:${image_id}"
		bootDiskSizeGb: 40
		cpu: bwa_cores
		memory: memory
		disks: "local-disk " + disk + " SSD"

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
		File mapped_pairs
		File mapped_stats
		String sample_id
		String memory = "20GB"
		String disk = "50"
		String image_id
	}

	command {
		python3 /home/qc_stats.py -i ${mapped_pairs} -p ${mapped_stats} -d ${sample_id} 
		zip -q ${sample_id}_qc.zip html_report.html hist.png
	}

	runtime {
		docker: "us-central1-docker.pkg.dev/aryeelab/docker/microc_qc:${image_id}"
		cpu: 1
		memory: memory
		disks: "local-disk " + disk + " SSD"
	}
	output {
		File qc_stats_file = "${sample_id}_qc.zip"
		String dist20kb = read_string(stdout())
	}
}
