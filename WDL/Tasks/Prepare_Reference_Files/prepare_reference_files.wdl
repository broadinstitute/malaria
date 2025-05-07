version 1.0

task prepare_reference_files {
	input {
		Array[String]? sample_ids
		File panel_info
		File reference_genome
		File? reference_amplicons
		File? reference_amplicons_2
		File? forward_primers_file
		File? reverse_primers_file
		File? path_to_snv

		Array[File] fastq1s
		Array[File] fastq2s
		File? barcodes_index
		File sample_metadata
	}
	command <<<
	export TMPDIR=tmp
	set -euxo pipefail
	
	###################################################################
	# Need to check for byte-order mark for interoperability purposes #
	###################################################################
	echo "Checking for byte-order mark (BOM) in panel_info file."
	has_bom() { head -c3 "$1" | grep -q $'\xef\xbb\xbf'; }
	if has_bom ~{panel_info}; then 
		echo "File has a BOM. Removing."
		sed -i 's/\xef\xbb\xbf//' ~{panel_info}
	else
		echo "File does not have a BOM. Skipping this step."
	fi

	echo "Finished checking for BOM."

	###################################################################
	# Make amplicon panel bed file from panel_info                    #
	###################################################################
	echo "Creating amplicon panel bed file."
	has_primers() { awk 'NR==1 {exit !(NF > 3)}' "$1"; }

	if [[ "~{forward_primers_file}" != '' && "~{reverse_primers_file}" != '' ]]; then
		echo "Custom primers provided. Using custom primers."
	elif has_primers ~{panel_info}; then
		echo "Panel info file contains primer information"
		cut -f1,2,3 ~{panel_info} | tail -n +2  | awk -v FS='\t' -v OFS='\t' '{print $1, $2-1, $3}' > amplicon_panel.bed	
	else
		echo "No primer files were provided. Please either provide primer information in panel_info or the forward_primers_file and reverse_primers_file."
		exit 1
	fi

	###################################################################
	# Make reference fasta file if reference not provided by user     #
	###################################################################
	echo "Checking for reference files."
	if [[ "~{reference_amplicons}" != '' ]]; then
		echo "Amplicon reference file provided."
		cp ~{reference_amplicons} reference.fasta
	elif [[ "~{reference_amplicons}" == '' && "~{reference_genome}" != '' && "~{panel_info}" != '' ]]; then
		echo "Amplicon reference file not provided. Reference genome and amplicon panel info file provided."
		echo "Creating amplicon reference file."
		bedtools getfasta -fi ~{reference_genome} -bed amplicon_panel.bed -fo reference_tmp.fasta
		echo "Created amplicon reference file."
		echo "Matching names to target_id in panel info."
		awk '
		BEGIN { FS=OFS="\t" }
		NR==FNR {
			if (FNR == 1) next;   # Skip header line in primers.tsv
			target_id[++i] = $6;
			next
		}
		/^>/ {
			print ">" target_id[++j];
			next
		}
		{ print }
		' "~{panel_info}" reference_tmp.fasta > reference.fasta	
	else
		echo "Neither reference file provided nor reference genome and amplicon panel info file provided."
		echo "Please provide any of these files to run the pipeline."
		exit 1
	fi

	echo "Finished checking for reference files."

	###################################################################
	# Make reference fasta file if reference not provided by user     #
	###################################################################
	echo "Checking for reference files."
	if [[ "~{reference_amplicons}" != '' ]]; then
		echo "Amplicon reference file provided."
		cp ~{reference_amplicons} reference.fasta
	elif [[ "~{reference_amplicons}" == '' && "~{reference_genome}" != '' && "~{panel_info}" != '' ]]; then
		echo "Amplicon reference file not provided. Reference genome and amplicon panel info file provided."
		echo "Creating amplicon reference file."
		bedtools getfasta -fi ~{reference_genome} -bed amplicon_panel.bed -fo reference.fasta
		echo "Created amplicon sequence file."
	else
		echo "Neither reference file provided nor reference genome and amplicon panel info file provided."
		echo "Please provide any of these files to run the pipeline."
	fi

	echo "Finished checking for reference files."

	###################################################################
	# Make forward and reverse primer files from panel_info           #
	###################################################################
	make_primers_from_panelinfo() {
		local input_tsv=$1
		local output_fasta=$2
		local seq_column=$3

		awk -v col="$seq_column" '
		BEGIN { FS=OFS="\t" }
		FNR == 1 { next }  # skip header
		{
		print ">" $6
		print $col
		}
		' "$input_tsv" > "$output_fasta"
	}

	echo "Checking for primer files."
	if [[ "~{forward_primers_file}" != '' ]]; then
		echo "Forward primers file provided."
		cp ~{forward_primers_file} primers_fw.fasta
	elif [[ "~{forward_primers_file}" == '' && "~{reference_genome}" != '' ]]; then
		if has_primers ~{panel_info}; then
			echo "Forward primers file not provided. Reference genome and amplicon panel info file provided."
			echo "Creating forward primers file."

			# Create forward primer fasta from panel_info
			make_primers_from_panelinfo ~{panel_info} primers_fw.fasta 4
			echo "Created forward primers file."
		else
			echo "Panel info file does not contain primer information"
			exit 1
		fi
	else
		echo "Neither forward primers file provided nor reference genome provided."
		echo "Please provide any of these files to run the pipeline."
	fi

	if [[ "~{reverse_primers_file}" != '' ]]; then
		echo "Reverse primers file provided."
		cp ~{reverse_primers_file} primers_rv.fasta
	elif [[ "~{reverse_primers_file}" == '' && "~{reference_genome}" != '' ]]; then
		if has_primers ~{panel_info}; then
			echo "Reverse primers file not provided. Reference genome and amplicon panel info file provided."
			echo "Creating reverse primers file."

			# Create reverse primer fasta from panel_info
			make_primers_from_panelinfo ~{panel_info} primers_rv.fasta 5
			echo "Created reverse primers file."
		else
			echo "Panel info file does not contain primer information"
			exit 1
		fi
	else
		echo "Neither reverse primers file provided nor reference genome and amplicon panel info file provided."
		echo "Please provide any of these files to run the pipeline."
	fi

	echo "Finished checking for primer files."

	#####################################################################
	# Enforce the same amplicon name across all files                   #
	#####################################################################
	# Before running in the next step, handle the special case where the 
	# files contain CI primers.
	# Function to extract the first 6 sequences with headers matching >KELT_
	extract_kelt_sequences() {
		local input_file=$1
		local output_file=$2
		local backup_file=$3

		# Backup original file
		cp "${input_file}" "${backup_file}"

		# Extract the first 6 sequences in the file (regardless of header)
		awk '
		BEGIN { RS=">"; ORS="" }
		NR > 1 { print ">" $0; count++ }
		count == 6 { exit }
		' "${input_file}" > "${output_file}"
	}

	# Run function in the forward and reverse primer files
	if grep -q ">KELT_" "~{forward_primers_file}"; then
		echo "CI tags detected"
		extract_kelt_sequences "~{forward_primers_file}" primers_fw.fasta primers_fw_CI.fasta
		extract_kelt_sequences "~{reverse_primers_file}" primers_rv.fasta primers_rv_CI.fasta
	fi

	echo "Checking for consistency across all panel reference files."
	python /Code/check_input_headers.py -i ~{panel_info} \
		--ref_genome ~{reference_genome} \
		~{if defined(reference_amplicons) then "--ref_amplicons_provided" else ""}

	#####################################################
	# Make bioinformatics_info.json			    #
	#####################################################

	fastq1s_string=$(IFS=" "; echo "~{sep=' ' fastq1s}")
	fastq1s_a=(${fastq1s_string})
	echo "ampseq.fastq1s" > samples_tmp.csv
	for i in "${!fastq1s_a[@]}"; do
		echo "${fastq1s_a[i]##*/}" >> samples_tmp.csv
	done

	python /Code/add_entry_to_json.py samples_tmp.csv config.json

	fastq2s_string=$(IFS=" "; echo "~{sep=' ' fastq2s}")
	fastq2s_a=(${fastq2s_string})
	echo "ampseq.fastq2s" > samples_tmp.csv
	for i in "${!fastq2s_a[@]}"; do
		echo "${fastq2s_a[i]##*/}" >> samples_tmp.csv
	done

	python /Code/add_entry_to_json.py samples_tmp.csv config.json

	echo "ampseq.panel_info" > samples_tmp.csv
	panel_info_tmp="~{panel_info}"
	echo "${panel_info_tmp##*/}" >> samples_tmp.csv
	python /Code/add_entry_to_json.py samples_tmp.csv config.json

	echo "ampseq.reference_genome" > samples_tmp.csv
	reference_genome_tmp="~{reference_genome}"
	echo "${reference_genome_tmp##*/}" >> samples_tmp.csv
	python /Code/add_entry_to_json.py samples_tmp.csv config.json

	echo "ampseq.sample_metadata" > samples_tmp.csv
	sample_metadata_tmp="~{sample_metadata}"
	echo "${sample_metadata_tmp##*/}" >> samples_tmp.csv
	python /Code/add_entry_to_json.py samples_tmp.csv config.json

	if [[ "~{reference_amplicons}" != '' ]]; then
		echo "ampseq.reference_amplicons" > samples_tmp.csv
		reference_amplicons_tmp="~{reference_amplicons}"
		echo "${reference_amplicons_tmp##*/}" >> samples_tmp.csv
		python /Code/add_entry_to_json.py samples_tmp.csv config.json
	fi

	if [[ "~{reference_amplicons_2}" != '' ]]; then
		echo "ampseq.reference_amplicons_2" > samples_tmp.csv
		reference_amplicons_2_tmp="~{reference_amplicons_2}"
		echo "${reference_amplicons_2_tmp##*/}" >> samples_tmp.csv
		python /Code/add_entry_to_json.py samples_tmp.csv config.json
	fi

	if [[ "~{forward_primers_file}" != '' ]]; then
		echo "ampseq.forward_primers_file" > samples_tmp.csv
		forward_primers_file_tmp="~{forward_primers_file}"
		echo "${forward_primers_file_tmp##*/}" >> samples_tmp.csv
		python /Code/add_entry_to_json.py samples_tmp.csv config.json
	fi

	if [[ "~{reverse_primers_file}" != '' ]]; then
		echo "ampseq.reverse_primers_file" > samples_tmp.csv
		reverse_primers_file_tmp="~{reverse_primers_file}"
		echo "${reverse_primers_file_tmp##*/}" >> samples_tmp.csv
		python /Code/add_entry_to_json.py samples_tmp.csv config.json
	fi

	if [[ "~{path_to_snv}" != '' ]]; then
		echo "ampseq.path_to_snv" > samples_tmp.csv
		path_to_snv_tmp="~{path_to_snv}"
		echo "${path_to_snv_tmp##*/}" >> samples_tmp.csv
		python /Code/add_entry_to_json.py samples_tmp.csv config.json
	fi

	if [[ "~{sep=' ' sample_ids}" != '' ]]; then
		echo "ampseq.sample_ids" > samples_tmp.csv
		sample_ids_string=$(IFS=" "; echo "~{sep=' ' sample_ids}")
		sample_ids_a=(${sample_ids_string})
		for i in "${!sample_ids_a[@]}"; do
			echo "${sample_ids_a[i]}" >> samples_tmp.csv
		done
		python /Code/add_entry_to_json.py samples_tmp.csv config.json
	fi

	if [[ "~{barcodes_index}" != '' ]]; then
		echo "ampseq.barcodes_index" > samples_tmp.csv
		barcodes_index_tmp="~{barcodes_index}"
		echo "${barcodes_index_tmp##*/}" >> samples_tmp.csv
		python /Code/add_entry_to_json.py samples_tmp.csv config.json
	fi

	>>>

	output {
		File reference_o = "reference.fasta"
		File panel_bedfile_o = "amplicon_panel.bed"
		File forward_primers_o = "primers_fw.fasta"
		File reverse_primers_o = "primers_rv.fasta"
		File? forward_primers_CI_o = "primers_fw_CI.fasta"
		File? reverse_primers_CI_o = "primers_rv_CI.fasta"
		File markers_table_o = "markers_table.csv"
		File bioinformatics_json_o = "config.json"
	}

	runtime {
		cpu: 1
		memory: "1 GB"
		disks: "local-disk 10 HDD"
		bootDiskSizeGb: 10
		preemptible: 3
		maxRetries: 1
		docker: "jorgeamaya/fileprep_ampseq:v_0_0_3"
	}
}

