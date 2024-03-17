version 1.0

workflow chunk_pe_fastq {

	input {
		String sample_id
		Array[File] fastq1
		Array[File] fastq2
		Int num_reads_per_chunk = 10000000
	}
	
	# Split the fastq files into chunks for parallelization
	scatter (fastq in fastq1) {
		call chunk_fastq_file as chunk_fastq1  { input:
					sample_id = sample_id,
					fastq = fastq,
					suffix = 'r1',
					num_lines_per_chunk = 4 * num_reads_per_chunk,
		}
	} 	

	scatter (fastq in fastq2) {
		call chunk_fastq_file as chunk_fastq2  { input:
					sample_id = sample_id,
					fastq = fastq,
					suffix = 'r2',
					num_lines_per_chunk = 4 * num_reads_per_chunk,
		}
	}
	
	call sort_chunks as sort_chunks1 { input:
		chunks = flatten(chunk_fastq1.chunks)
	}

	call sort_chunks as sort_chunks2 { input:
		chunks = flatten(chunk_fastq2.chunks)
	}
	
	output {
		Array[File] fastq1_chunks = sort_chunks1.sorted_chunks
		Array[File] fastq2_chunks = sort_chunks2.sorted_chunks
	}
}

task chunk_fastq_file {
    input {
		String sample_id
        File fastq
        String suffix
        Int num_lines_per_chunk
        Int disk_gb = 375
        Int disk_gb_local = ceil(disk_gb/375) * 375
    }
    
    command {
        pigz -dc -p2 ${fastq} | split -d --suffix-length=3 -l ${num_lines_per_chunk} --additional-suffix='_${suffix}.fastq' --filter='pigz -c -p12 > $FILE.gz' - ${sample_id}-
    }
    
     runtime {
        continueOnReturnCode: false
        docker: "us-central1-docker.pkg.dev/aryeelab/docker/utils"
        cpu: 16
        disks: "local-disk " + disk_gb_local + " LOCAL"        
    }   
    output {
        Array[File] chunks = glob("*.fastq.gz")
    }
}

task sort_chunks { 
	input {
		Array[String] chunks
	}
	
	command {
		echo ${sep=" " chunks} | tr " " "\n" | sort
	}
	
	runtime {
        continueOnReturnCode: false
        docker: "ubuntu"
        cpu: 1
    } 
    
    output {
        Array[String] sorted_chunks = read_lines(stdout())
    }
    
}