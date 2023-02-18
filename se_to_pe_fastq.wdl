version 1.0

workflow se_to_pe_fastq {

	input {
		Array[File] fastq
		Int num_lines_per_chunk = 40000000
	}

	scatter (fq in fastq) {
		call se_to_pe {input: fastq = fq}
	}
	
	output {
		Array[File] fastq_r1 = se_to_pe.fastq_r1
		Array[File] fastq_r2 = se_to_pe.fastq_r2
	}
}

task se_to_pe {
	input {
		File fastq
		String read1_length = 200
		String read2_length = 100		
		String fq1 = basename(fastq, ".fastq.gz") + ".r1.fastq.gz"
		String fq2 = basename(fastq, ".fastq.gz") + ".r2.fastq.gz"
	}

	command {		
		seqkit subseq -r 1:${read1_length} ${fastq} -o ${fq1}
		seqkit seq -r -p ${fastq} | seqkit subseq -r 1:${read2_length} -o ${fq2} 
	}
	
	runtime {
		docker: "us-central1-docker.pkg.dev/aryeelab/docker/seqkit"
		cpu: 8
		memory: "8G"
		disks: "local-disk " + 500 + " SSD" 
	}

	output {
		File fastq_r1 = "${fq1}"
		File fastq_r2 = "${fq2}"
	}
}




