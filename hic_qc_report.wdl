workflow hic_qc_report {
    Array[File] qc_stats_zips
    call read_stats_report {input: qc_stats_zips = qc_stats_zips}
}


task read_stats_report {
        Array[File] qc_stats_zips

        command <<<
                            
            mkdir hicpro_stats
            cd hicpro_stats
            for z in ${sep=' ' qc_stats_zips}; do
                unzip $z
            done
            cd ..
            
            Rscript /usr/local/bin/HiCStats.R -d hicpro_stats -o qc

        >>>
                
        output {
            File cisTrans = "qc-cisTrans.png"
            File numberOfPairs = "qc-numberOfPairs.png"
            File pairValidity = "qc-pairValidity.png"
        }
        
        
        runtime {
            continueOnReturnCode: false
            docker: "aryeelab/hic_qc_report:latest"
        }
}
