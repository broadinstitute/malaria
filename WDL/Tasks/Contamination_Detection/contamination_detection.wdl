version 1.0

task contamination_detection {
	input {
		Array[File] adaptor_rem1s
		Array[File] adaptor_rem2s
		File? forward_primers_file
		File? reverse_primers_file
		File? barcodes_index

		#Contamination detection parameters
		Int minreads_threshold = 1000
		Float contamination_threshold = 1
	}

	command <<<
	export TMPDIR=tmp
	set -euxo pipefail

	###################################################################
	# Contamination detection					  #
	###################################################################
	echo "Performing contamination detection"
	
	adaptor_rem1s_string=$(IFS=" "; echo "~{sep=' ' adaptor_rem1s}")
	adaptor_rem2s_string=$(IFS=" "; echo "~{sep=' ' adaptor_rem2s}")
	
	# Detect missing files
	touch missing_files.tsv
	touch found_files.tsv

	cut -d, -f1 "~{barcodes_index}" | while read -r sample; do
		if [ "$sample" == "sample_id" ]; then
			continue
		fi

		if ! echo "$adaptor_rem1s_string" | grep -qi "$sample"; then
			echo "$sample" >> missing_files.tsv
		else
			echo "$sample" >> found_files.tsv
		fi
		echo "$sample" >> all_files.tsv
	done

	echo "Sequencing run with inline barcodes. Performing analysis of combinatorial indices."
	python /Code/CI_TerraPipeline.py -f1 "${adaptor_rem1s_string}" -f2 "${adaptor_rem2s_string}" -b "~{barcodes_index}" -l found_files.tsv

	Rscript /Code/render_report.R -d $PWD/Report/Merge/ -o $PWD/Report/ -p ~{barcodes_index} -m 1000 -c 0.5 -mf $PWD/missing_files.tsv
	tar -czvf Report_Cards.tar.gz Report
	
	echo "Finished contamination detection"

>>>

	output {
		File missing_files = "missing_files.tsv"
		File? contamination_detection_sample_cards = "Report_Cards.tar.gz"
		File? contamination_detection_report = "Results/ci_report_layouting.html"
	}
	runtime {
		cpu: 1
		memory: "15 GiB"
		disks: "local-disk 10 HDD"
		bootDiskSizeGb: 10
		preemptible: 3
		maxRetries: 1
		docker: 'jorgeamaya/ci_processing:v_0_0_2'
	}
}
