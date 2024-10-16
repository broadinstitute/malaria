version 1.0

task trimprimers {
	input {
		File trimmed_fastq1
		File trimmed_fastq2
		File forward_primers_file
		File reverse_primers_file
	}

	command <<<
	export TMPDIR=tmp
	set -euxo pipefail

	echo "Removing primers"
	mkdir Results

	sample=$(basename "~{trimmed_fastq1}" | perl -pe "s/_val_1\.fq\.gz//")

	cutadapt -g file:~{forward_primers_file} \
		-G file:~{reverse_primers_file} \
		-o Results/${sample}_mixed_op_1.fq.gz \
		-p Results/${sample}_mixed_op_2.fq.gz \
		--pair-adapters --action=trim \
		--discard-untrimmed \
		--cores=5 \
		~{trimmed_fastq1} ~{trimmed_fastq2}

	>>>
	output {
		Array[File] primerrem = glob("Results/*")
	}
	runtime {
		cpu: 5
		memory: "5 GiB"
		disks: "local-disk 10 HDD"
		bootDiskSizeGb: 10
		preemptible: 0
		maxRetries: 1
		docker: 'jorgeamaya/trimprimers:testing'
	}
}
