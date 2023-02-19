version 1.0

workflow download_url_to_file {

	input {
		Array[String] urls
	}
	
	scatter (url in urls) {
		call download {input: url = url}
	}
	
	
	output {
		Array[File] file = download.file
	}

}

task download {
	input {
		String url
		String filename = basename(url)
	}

	command {
		echo "Downloading ${url}..."
		curl -o ${filename} ${url}
	}
	
	runtime {
		docker: "us-central1-docker.pkg.dev/aryeelab/docker/utils"
		cpu: 2
		memory: "4G"
		disks: "local-disk " + 500 + " SSD" 
	}

	output {
		File file = "${filename}"
	}
}




