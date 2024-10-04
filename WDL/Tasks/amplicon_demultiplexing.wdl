version 1.0

#task ampseq_pipeline_demult {
#	input {
#		Array[File] path_to_r1
#		Array[File] path_to_r2
#		File path_to_flist
#		File pr1
#		File pr2
#		File reference
#		File? reference2
#		File? path_to_snv
#		File config_json
#	}
#
#	command <<<
#	set -euxo pipefail
#	cat ~{config_json}
#
#	./Code/file_staging.sh
#
#	echo "Demultiplexing requested."
#	echo "Demultiplexing will be performed in the data. Some read pairs assumed to be too short to overlap and be merged."
#
#	#Check if the first line in barcodes_matches.csv indicates the presence of inline barcodes
#	if grep -q "," references/barcodes_matches ; then
#		echo "Sequencing run with inline barcodes. Performing analysis of combinatorial indices followed by denoising"
#		find . -type f
#		python /Code/Amplicon_TerraPipeline.py --config ~{config_json} --terra --meta --adaptor_removal --contamination --demultiplex --primer_removal
#		find . -type f
#		Rscript /Code/render_report.R -d /cromwell_root/Report/Merge/ -o /cromwell_root/Report/ -p ~{path_to_flist} -m 1000 -c 0.5 -mf /cromwell_root/Results/missing_files.tsv
#		find . -type f
#		tar -czvf Report_Cards.tar.gz /cromwell_root/Report
#		find . -type f
#	else
#		echo "Sequencing run without inline barcodes. Skipping analysis of combinatorial indices and performing only denoising"
#		find . -type f
#		python /Code/Amplicon_TerraPipeline.py --config ~{config_json} --terra --meta --adaptor_removal --demultiplex --primer_removal
#		find . -type f
#	fi
#	
#	>>>
#	output {
#		File missing_files = "Results/missing_files.tsv"
#		Array[File] PrimerRem = glob("Results/PrimerRem/*")
#		Array[File] AdaptorRem = glob("Results/AdaptorRem/*")
#		File? decontamination_sample_cards = "Report_Cards.tar.gz"
#		File? decontamination_report = "Results/ci_report_layouting.html"
#	}
#	runtime {
#		cpu: 1
#		memory: "15 GiB"
#		disks: "local-disk 10 HDD"
#		bootDiskSizeGb: 10
#		preemptible: 3
#		maxRetries: 1
#		docker: 'jorgeamaya/ampseq:latest'
#	}
#}
#
