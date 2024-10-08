version 1.0

task amplicon_no_demultiplexing {
	input {
		Array[File] fastq1s
		Array[File] fastq2s
		File path_to_flist
		File forward_primers_file
		File reverse_primers_file
		File reference
		File? reference2
		File? path_to_snv
		File config_json
	}

	command <<<
	set -euxo pipefail
	cat ~{config_json}

	###################################################################
	# Copy files to appropriate folders within Docker                 #
	###################################################################

	mkdir fq_dir
	mkdir references
	gsutil -m cp -r ~{sep = ' ' fastq1s} fq_dir/
	gsutil -m cp -r ~{sep = ' ' fastq2s} fq_dir/
	gsutil cp ~{path_to_flist} references/samples.txt
	gsutil -m cp -r ~{forward_primers_file} references/
	gsutil -m cp -r ~{reverse_primers_file} references/
	gsutil -m cp -r ~{reference} references/

	~{"gsutil cp " + reference2 + " references/reference2.fasta"}
	~{"gsutil cp " + path_to_snv + " references/snv_filter.tsv"}

	echo "Demultiplexing not requested."
	echo "No demultiplexing will be performed in the data. Read pairs assumed to be long enough to overlap and be merged."

	echo "Sequencing run without inline barcodes. Skipping analysis of combinatorial indices."
	find . -type f
	python /Code/Amplicon_TerraPipeline.py --config ~{config_json} --terra --meta --adaptor_removal --primer_removal
	find . -type f
	>>>
	output {
		File missing_files = "Results/missing_files.tsv"
		Array[File] PrimerRem = glob("Results/PrimerRem/*")
		Array[File] AdaptorRem = glob("Results/AdaptorRem/*")
		File? decontamination_sample_cards = "Report_Cards.tar.gz"
		File? decontamination_report = "Results/ci_report_layouting.html"
	}
	runtime {
		cpu: 1
		memory: "15 GiB"
		disks: "local-disk 10 HDD"
		bootDiskSizeGb: 10
		preemptible: 0
		maxRetries: 1
		docker: 'jorgeamaya/ampseq:latest'
	}
}
