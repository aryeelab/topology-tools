version 1.0

workflow se_to_pe_fastq {

	input {
		File fastq
	}

	call se_to_pe {input: fastq = fastq}
	
	output {
		File fastq1 = se_to_pe.fastq1
		File fastq2 = se_to_pe.fastq2
	}

}

task se_to_pe {
	input {
		File fastq
		String read_length = 150
		String fq1 = basename(fastq, ".fastq.gz") + ".r1.fastq.gz"
		String fq2 = basename(fastq, ".fastq.gz") + ".r2.fastq.gz"
	}

	command {		
		seqkit subseq -r 1:${read_length} ${fastq} | gzip -c  > ${fq1}
		seqkit seq -r -p ${fastq} | seqkit subseq -r 1:${read_length} | gzip -c  > ${fq2} 
	}
	
	runtime {
		docker: "us-central1-docker.pkg.dev/aryeelab/docker/seqkit"
		cpu: 8
		memory: "8G"
		disks: "local-disk " + 2000 + " SSD" 
	}

	output {
		File fastq1 = "${fq1}"
		File fastq2 = "${fq2}"
	}
}




