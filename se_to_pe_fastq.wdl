version 1.0

workflow se_to_pe_fastq {

	input {
		File fastq
	}

  # Split the fastq file into chunks for parallelization
	call chunk_fastq_file  { input: 
		fastq = fastq,
		num_lines_per_chunk = 40000000
	}

	scatter (fastq in chunk_fastq_file.chunks) {
		call se_to_pe {input: fastq = fastq}
	}
	
	output {
		Array[Pair[File, File]] fastq_pairs = se_to_pe.fastq_pairs
	}

}

task chunk_fastq_file {
    input {
        File fastq
        String sample_id = basename(fastq, ".fastq.gz")
        Int num_lines_per_chunk
    }
    
    command {
        pigz -dc -p2 ${fastq} | split -d --suffix-length=3 -l ${num_lines_per_chunk} --filter='pigz -c -p24 > $FILE.fastq.gz' - ${sample_id}-
    }
    
     runtime {
        continueOnReturnCode: false
        docker: "us-central1-docker.pkg.dev/aryeelab/docker/utils"
        cpu: 32
        disks: "local-disk 375 LOCAL"        
    }   
    output {
        Array[File] chunks = glob("*.fastq.gz")
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
		seqkit subseq -r 1:${read_length} ${fastq} -o ${fq1}
		seqkit seq -r -p ${fastq} | seqkit subseq -r 1:${read_length} -o ${fq2} 
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
		Pair[File, File] fastq_pairs = (fastq1, fastq2)
	}
}




