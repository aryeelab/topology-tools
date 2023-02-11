version 1.0

workflow download_url_to_file {

	call download {}
	
	output {
		File file = download.file
	}

}

task download {
	input {
		String basename
		String extension = "fastq.gz"
		String urls
		String out_file = basename + '.' + extension
		String arr = "{ADDR[@]}" 
	}

	command {
		IFS=',' read -ra ADDR <<< "${urls}"
		for URL in "$${arr}"; do
			echo "Processing $URL"
			curl $URL >> ${out_file}
		done	
	}
	
	runtime {
		docker: "ubuntu"
		cpu: 2
		memory: "4G"
		disks: "local-disk " + 300 + " SSD" 
	}


	output {
		File file = "${out_file}"
	}
}




