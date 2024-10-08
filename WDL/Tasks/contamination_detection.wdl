version 1.0

task contamination_detection {
	input {
		Array[File] fastq1s
		Array[File] fastq2s
		File path_to_flist
		File forward_primers_file
		File reverse_primers_file
		File config_json
		File? barcodes_index
	}

	command <<<
	set -euxo pipefail
	cat ~{config_json}

	mkdir fq_dir
	mkdir references
	gsutil -m cp -r ~{sep = ' ' fastq1s} fq_dir/
	gsutil -m cp -r ~{sep = ' ' fastq2s} fq_dir/
	gsutil cp ~{path_to_flist} references/samples.txt
	gsutil -m cp -r ~{forward_primers_file} references/
	gsutil -m cp -r ~{reverse_primers_file} references/
	gsutil cp ~{barcodes_index} references/barcodes_index.csv

	echo "Performing contamination detection"

	###################################################################
	# Check if barcodes are present in barcodes.csv file              #
	###################################################################

	echo "Sequencing run with inline barcodes. Performing analysis of combinatorial indices."
	find . -type f
	python /Code/CI_TerraPipeline.py --config ~{config_json} --terra --meta --adaptor_removal --contamination
	find . -type f
	Rscript /Code/render_report.R -d /cromwell_root/Report/Merge/ -o /cromwell_root/Report/ -p ~{path_to_flist} -m 1000 -c 0.5 -mf /cromwell_root/Results/missing_files.tsv
	find . -type f
	tar -czvf Report_Cards.tar.gz /cromwell_root/Report
	find . -type f
	>>>
	output {
		File missing_files = "Results/missing_files.tsv"
		File? contamination_detection_sample_cards = "Report_Cards.tar.gz"
		File? contamination_detection_report = "Results/ci_report_layouting.html"
	}
	runtime {
		cpu: 1
		memory: "15 GiB"
		disks: "local-disk 10 HDD"
		bootDiskSizeGb: 10
		preemptible: 0
		maxRetries: 1
		docker: 'jorgeamaya/ci_processing:latest'
	}
}
