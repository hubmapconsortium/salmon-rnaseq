version 1.0

task fastqc {
	input {
		Array[File] fastqs
		File? limits
		Int addldisk = 10
		Int threads = 4
		Int mem_gb = 8
		Int preempt = 1 
	}
	Int finalDiskSize = addldisk + ceil(size(fastqs, "GB"))

	command {
		mkdir outputs
		if [[ "~{limits}" != "" ]]
		then
			fastqc -o outputs -l ~{limits} ~{sep=" " fastqs}
		else
			fastqc -o outputs ~{sep=" " fastqs}
		fi
	}

	runtime {
		#threads: threads
		docker: "hubmap/scrna-analysis:latest"
		#disks: "local-disk " + finalDiskSize + " SSD"
		#memory: "${mem_gb} GB"
		#preemptible: "${preempt}"
	}

	output {
		Array[File] reports = glob("outputs/*.html")
	}	
}