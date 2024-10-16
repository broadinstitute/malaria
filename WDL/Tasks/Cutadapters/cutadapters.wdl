version 1.0

task cutadapters {
	input {
		File fastq1
		File fastq2
		Int trim_galore_qvalue
		Int trim_galore_length
		String pattern_fw
		String pattern_rv
	}

	command <<<
	export TMPDIR=tmp
	set -euxo pipefail
	
	echo "Removing adapters"
	mkdir Results

	sample=$(basename "~{fastq1}" | perl -pe "s/~{pattern_fw}//")
	
	trim_galore --paired --gzip \
		--quality ~{trim_galore_qvalue} \
		--length ~{trim_galore_length} \
		--output_dir Results \
		--basename ${sample} \
		~{fastq1} ~{fastq2}
	>>>

	output {
		# Generate a list of trimmed forward and reverse files
		Array[File] fastq1s_noadapters = glob("Results/*_val_1.fq.gz")
		Array[File] fastq2s_noadapters = glob("Results/*_val_2.fq.gz") 
	}

	runtime {
		cpu: 1
		memory: "5 GiB"
		disks: "local-disk 10 HDD"
		bootDiskSizeGb: 10
		preemptible: 3
		maxRetries: 1
		docker: 'jorgeamaya/cutadapters:testing'
	}
}
