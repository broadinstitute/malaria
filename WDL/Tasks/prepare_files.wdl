version 1.0

task prepare_files {
	input {
		Array[String] sample_ids
		File amplicon_info

		File? reference_amplicons
		File? reference_amplicons_2

		File? forward_primers_file
		File? reverse_primers_file

		File? path_to_snv
		File? reference_genome
		String pattern_fw = "*_L001_R1_001.fastq.gz"
		String pattern_rv = "*_L001_R2_001.fastq.gz"
		String Class = "parasite"
		String maxEE = "5,5"
		String trimRight = "0,0"
		Int minLen = 30
		String truncQ = "5,5"
		String matchIDs = "0"
		Int max_consist = 10
		Float omegaA = 0.000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001
		String saveRdata = ""
		Int justConcatenate = 0
		Int maxMismatch = 0
		String no_ref = 'False'
		String adjust_mode = "absolute"
		String strain = "3D7"
		String strain2 = "DD2"
		String polyN = "5"
		String min_reads = "0"
		String min_samples = "0"
		String max_snv_dist = "-1"
		String max_indel_dist = "-1"
		String include_failed = "False"
		String exclude_bimeras = "False"
		Int minreads_threshold = 1000
		Float contamination_threshold = 0.5
		String verbose = "False"
		String adapter = "None"
	}
	
	File path_to_flist = write_lines(sample_ids)
	
	String forward_primers_file_refdir = "primers_fw.fasta" 
	String reverse_primers_file_refdir = "primers_rv.fasta"

	Map[String, String] in_map = {
		"path_to_fq": "fq_dir",
		"path_to_flist": "references/samples.txt",
		"forward_primers_file": forward_primers_file_refdir,
		"reverse_primers_file": reverse_primers_file_refdir,
		"reference_amplicons": "references/reference.fasta",
		"pattern_fw": pattern_fw,
		"pattern_rv": pattern_rv,
		"Class": Class,
		"maxEE": maxEE,
		"trimRight": trimRight,
		"minLen": minLen,
		"truncQ": truncQ,
		"matchIDs": matchIDs,
		"max_consist": max_consist,
		"omegaA": omegaA,
		"saveRdata": saveRdata,
		"justConcatenate": justConcatenate,
		"maxMismatch": maxMismatch,
		"no_ref": no_ref,
		"adjust_mode": adjust_mode,
		"strain": strain,
		"strain2": strain2,
		"polyN": polyN,
		"min_reads": min_reads,
		"min_samples": min_samples,
		"max_snv_dist": max_snv_dist,
		"max_indel_dist": max_indel_dist,
		"include_failed": include_failed,
		"exclude_bimeras": exclude_bimeras,
		"minreads_threshold": minreads_threshold,
		"contamination_threshold": contamination_threshold,
		"verbose": verbose,
		"adapter": adapter
	}

	File config_json = write_json(in_map)

	command <<<
	set -euxo pipefail

	###################################################################
	# Need to check for byte-order mark for interoperability purposes #
	###################################################################

	has_bom() { head -c3 "$1" | grep -q $'\xef\xbb\xbf'; }
	if has_bom ~{amplicon_info}; then 
		echo "File has a BOM. Removing."
		sed -i 's/\xef\xbb\xbf//' ~{amplicon_info}
	else
		echo "File does not have a BOM. Skipping this step."
	fi

	mv ~{path_to_flist} samples.txt

	###################################################################
	# Make reference fasta file if reference not provided by user     #
	###################################################################
	cut -f1,4,5 ~{amplicon_info} | tail -n +2 > amplicon_panel.bed

	if [[ "~{reference_amplicons}" != '' ]]; then
		echo "Amplicon reference file provided."
		cp ~{reference_amplicons} reference.fasta
	elif [[ "~{reference_amplicons}" == '' && "~{reference_genome}" != '' && "~{amplicon_info}" != '' ]]; then
		echo "Amplicon reference file not provided. Reference genome and amplicon panel info file provided."
		echo "Creating amplicon reference file."
		bedtools getfasta -fi ~{reference_genome} -bed amplicon_panel.bed -fo reference.fasta
		echo "Created amplicon sequence file."
	else
		echo "Neither reference file provided nor reference genome and amplicon panel info file provided."
		echo "Please provide any of these files to run the pipeline."
	fi

	###################################################################
	# Make forward and reverse primer files from amplicon_info        #
	###################################################################

	if [[ "~{forward_primers_file}" != '' ]]; then
		echo "Forward primers file provided."
		cp ~{forward_primers_file} primers_fw.fasta
	elif [[ "~{forward_primers_file}" == '' && "~{reference_genome}" != '' && "~{amplicon_info}" != '' ]]; then
		echo "Forward primers file not provided. Reference genome and amplicon panel info file provided."
		echo "Creating forward primers file."

		# Create forward primer fasta from amplicon_info
		awk -F'\t' '{print $1 "\t" $2 "\t" $4-1}' ~{amplicon_info} | tail -n +2 > primers_fw.bed
		bedtools getfasta -fo primers_fw.fasta -fi ~{reference_genome} -bed primers_fw.bed
		
		echo "Created forward primers file."
		rm primers_fw.bed
	else
		echo "Neither forward primers file provided nor reference genome and amplicon panel info file provided."
		echo "Please provide any of these files to run the pipeline."
	fi

	if [[ "~{reverse_primers_file}" != '' ]]; then
		echo "Reverse primers file provided."
		cp ~{reverse_primers_file} primers_rv.fasta
	elif [[ "~{reverse_primers_file}" == '' && "~{reference_genome}" != '' && "~{amplicon_info}" != '' ]]; then
		echo "Reverse primers file not provided. Reference genome and amplicon panel info file provided."
		echo "Creating reverse primers file."

		# Create reverse primer fasta from amplicon_info
		sed 's/\r//' ~{amplicon_info} | awk -F"\t" '{print $1 "\t" $5 "\t" $3 "\t\t\t-"}' | tail -n +2 > primers_rv.bed
		bedtools getfasta -fo primers_rv.fasta -fi ~{reference_genome} -bed primers_rv.bed
		
		echo "Created reverse primers file."
		rm primers_rv.bed
	else
		echo "Neither reverse primers file provided nor reference genome and amplicon panel info file provided."
		echo "Please provide any of these files to run the pipeline."
	fi



	###################################################################
	# Add reference_amplicons_2 and path_to_snv into config_json      #
	###################################################################
	if [[ "~{reference_amplicons_2}" != '' ]]; then
		echo "Reference amplicon 2 provided"
		mv ~{reference_amplicons_2} reference2.fasta
		python /Code/add_entry_to_json.py ~{config_json} "reference_amplicons_2" "reference2.fasta"
	fi

	if [[ "~{path_to_snv}" != '' ]]; then
		echo "Path to SNV file provided"
		mv ~{path_to_snv} snv_filter.tsv
		python /Code/add_entry_to_json.py ~{config_json} "path_to_snv" "snv_filter.tsv"
	fi
	cat ~{config_json} >> config_json.json

	>>>

	output {
		File config_json_out = "config_json.json"
		File path_to_flist_o = "samples.txt"
		File reference_out = "reference.fasta"
		File panel_bedfile = "amplicon_panel.bed"

		File? forward_primers_o = "primers_fw.fasta"
		File? reverse_primers_o = "primers_rv.fasta"

	}

	runtime {
		cpu: 1
		memory: "1 GB"
		disks: "local-disk 10 HDD"
		bootDiskSizeGb: 10
		preemptible: 3
		maxRetries: 1
		docker: "jorgeamaya/fileprep_ampseq:latest"
	}
}
