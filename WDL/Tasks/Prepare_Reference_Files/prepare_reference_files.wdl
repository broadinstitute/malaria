version 1.0

task prepare_reference_files {
	input {
		File panel_info
		File reference_genome
		File? reference_amplicons
		File? reference_amplicons_2
		File? forward_primers_file
		File? reverse_primers_file
		File? path_to_snv
	}
	command <<<
	# FUTURE DEVELOPMENT: 
	# 1 MAKE REFERENCE_OUT OPTIONAL. SKIP THE GENERATION OF THE REFERENCE FILE IF PROVIDED BY USER.
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
	if has_primers ~{panel_info}; then
		echo "Panel info file contains primer information"
		cut -f1,4,5 ~{panel_info} | tail -n +2  | awk -v FS='\t' -v OFS='\t' '{print $1, $2-1, $3}' > amplicon_panel.bed
	else
		echo "Panel info file does not contain primer information."
		if [[ "~{forward_primers_file}" == ''  || "~{reverse_primers_file}" == '' ]]; then
			echo "No primer files were provided. Please either provide primer information in panel_info or the forward_primers_file and reverse_primers_file."
			exit 1
		else
			tail -n +2 ~{panel_info} | awk -v FS='\t' -v OFS='\t' '{print $1, $2-1, $3}' > amplicon_panel.bed	
		fi
	fi

	echo "Created amplicon panel bed file."

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
	echo "Checking for primer files."
	if [[ "~{forward_primers_file}" != '' ]]; then
		echo "Forward primers file provided."
		cp ~{forward_primers_file} primers_fw.fasta
	elif [[ "~{forward_primers_file}" == '' && "~{reference_genome}" != '' ]]; then
		if has_primers ~{panel_info}; then
			echo "Forward primers file not provided. Reference genome and amplicon panel info file provided."
			echo "Creating forward primers file."

			# Create forward primer fasta from panel_info
			awk -F'\t' '{print $1 "\t" $2 "\t" $4-1}' ~{panel_info} | tail -n +2 > primers_fw.bed
			bedtools getfasta -fo primers_fw.fasta -fi ~{reference_genome} -bed primers_fw.bed
			
			echo "Created forward primers file."
			rm primers_fw.bed
		else
			echo "Panel info file does not contain primer information"
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
			sed 's/\r//' ~{panel_info} | awk -F"\t" '{print $1 "\t" $5 "\t" $3 "\t\t\t-"}' | tail -n +2 > primers_rv.bed
			bedtools getfasta -s -fo primers_rv.fasta -fi ~{reference_genome} -bed primers_rv.bed #Reverse complement

			echo "Created reverse primers file."
			rm primers_rv.bed
		else
			echo "Panel info file does not contain primer information"
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
	>>>

	output {
		File reference_o = "reference.fasta"
		File panel_bedfile_o = "amplicon_panel.bed"
		File forward_primers_o = "primers_fw.fasta"
		File reverse_primers_o = "primers_rv.fasta"
		File? forward_primers_CI_o = "primers_fw_CI.fasta"
		File? reverse_primers_CI_o = "primers_rv_CI.fasta"
		File markers_table_o = "markers_table.csv"
	}

	runtime {
		cpu: 1
		memory: "1 GB"
		disks: "local-disk 10 HDD"
		bootDiskSizeGb: 10
		preemptible: 3
		maxRetries: 1
		docker: "jorgeamaya/fileprep_ampseq:v_0_0_2"
	}
}

