version 1.0


workflow RCMC {
	String pipeline_ver = 'dev'
	String image_id = sub(pipeline_ver, "dev", "latest")

	input {
		File pairs
		File chrom_sizes
		String region
		String sample_id
		String hic_resolution
	}

	call filter_mapped_pairs_region { input:
		pairs = pairs,
        region = region,
        sample_id = sample_id
	}

	call juicer_hic {input: 
						image_id = image_id, 
						sample_id = sample_id, 
						chrom_sizes = chrom_sizes, 
						mapped_pairs = filter_mapped_pairs_region.captured_pairs
					}

	call juicer_hic_res {input: 
						image_id = image_id, 
						sample_id = sample_id, 
						chrom_sizes = chrom_sizes, 
						mapped_pairs = filter_mapped_pairs_region.captured_pairs,
						hic_resolution = hic_resolution
					}

	call cooler_res {input: 
						image_id = image_id, 
						sample_id = sample_id, 
						chrom_sizes = chrom_sizes, 
						mapped_pairs = filter_mapped_pairs_region.captured_pairs,
						hic_resolution = hic_resolution
					}

	output {
		File captured_pairs = filter_mapped_pairs_region.captured_pairs
		String captured_perc = filter_mapped_pairs_region.perc
		File hic = juicer_hic.hic
		File hic_res = juicer_hic_res.hic
		File cooler_raw_res = cooler_res.raw_mcool
		File cooler_balanced_res = cooler_res.balanced_mcool		
	}

}

task filter_mapped_pairs_region {
    input {
        File pairs
        String region
        String sample_id
        String disk = "20"
        String memory = "40GB"
    }
    command {
        python /home/RCMC.py -i ${pairs} -r ${region} -o ${sample_id}_captured.pairs
    }
    runtime {
        docker: "salvacasani/rcmc:latest"
        bootDiskSizeGb: 40
		memory: memory
		disks: "local-disk " + disk + " SSD"
    }
    output {
        String perc = read_string(stdout())
        File captured_pairs = "${sample_id}_captured.pairs"
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
			${sample_id}_captured.hic \
			${chrom_sizes}
	}

	runtime {
		docker: "us-central1-docker.pkg.dev/aryeelab/docker/juicer:${image_id}"
		bootDiskSizeGb: 40
		memory: memory
		disks: "local-disk " + disk + " SSD"
	}

	output {
		File hic = "${sample_id}_captured.hic"
	
	}

}


task juicer_hic_res {
	input {
		String image_id
		String sample_id
		File chrom_sizes
		File mapped_pairs
		Int cores = 2
		String memory = "40GB"
		String disk = "200"
		String hic_resolution
	}

	command {
		java -Xmx32000m  -Djava.awt.headless=true -jar /usr/local/bin/juicer_tools_1.22.01.jar pre \
			--threads ${cores} -r ${hic_resolution}\
			${mapped_pairs} \
			${sample_id}_${hic_resolution}bp_captured.hic \
			${chrom_sizes}
	}

	runtime {
		docker: "us-central1-docker.pkg.dev/aryeelab/docker/juicer:${image_id}"
		bootDiskSizeGb: 40
		memory: memory
		disks: "local-disk " + disk + " SSD"
	}

	output {
		File hic = "${sample_id}_${hic_resolution}bp_captured.hic"
	
	}

}

task cooler_res {
	input {
		String image_id
		String sample_id
		File chrom_sizes
		File mapped_pairs
		Int hic_resolution = "500"
		String disk = "200"
	}

	command {
		cooler cload pairs -c1 2 -p1 3 -c2 4 -p2 5 ${chrom_sizes}:${hic_resolution} ${mapped_pairs} ${sample_id}.cool
		cooler zoomify --resolutions ${hic_resolution}N -o ${sample_id}.captured.raw.mcool -p 4 ${sample_id}.cool
		cooler zoomify --resolutions ${hic_resolution}N -o ${sample_id}.captured.balanced.mcool -p 4 --balance --balance-args '--nproc 4' ${sample_id}.cool
	}

	runtime {
		docker: "us-central1-docker.pkg.dev/aryeelab/docker/cooler:${image_id}"
		cpu: 4
		memory: "8GB"
		disks: "local-disk " + disk + " SSD"
	}

	output {
		File raw_mcool = "${sample_id}.captured.raw.mcool"	
		File balanced_mcool = "${sample_id}.captured.balanced.mcool"			
	}
}

