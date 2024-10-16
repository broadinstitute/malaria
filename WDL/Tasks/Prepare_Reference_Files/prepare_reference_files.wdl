version 1.0

task prepare_reference_files {
	input {
		Array[String] sample_ids
		File forward_primers_file
		File reverse_primers_file
		File? reference_amplicons
		File? reference_amplicons_2
		File? path_to_snv
		File? reference_genome
		File? panel_bedfile
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
		File? barcodes_index
	}
	File path_to_flist = write_lines(sample_ids)
	String forward_primers_file_refdir = "references/" + basename(forward_primers_file)
	String reverse_primers_file_refdir = "references/" + basename(reverse_primers_file)

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
	export TMPDIR=tmp
	set -euxo pipefail

	###################################################################
	# Need to check for byte-order mark for interoperability purposes #
	###################################################################

	has_bom() { head -c3 "$1" | grep -q $'\xef\xbb\xbf'; }
	if has_bom ~{panel_bedfile}; then 
		echo "File has a BOM. Removing."
		sed -i 's/\xef\xbb\xbf//' ~{panel_bedfile}
	else
		echo "File does not have a BOM. Skipping this step."
	fi

	###################################################################
	# Make reference fasta file if reference not provided by user     #
	###################################################################

	if [[ "~{reference_amplicons}" != '' ]]; then
		echo "Reference file provided."
		cp ~{reference_amplicons} reference.fasta
	elif [[ "~{reference_amplicons}" == '' && "~{reference_genome}" != '' && "~{panel_bedfile}" != '' ]]; then
		echo "Reference file not provided."
		echo "Reference genome and panel's bed file provided."
		echo "Creating reference file"
		bedtools getfasta -fi ~{reference_genome} -bed ~{panel_bedfile} -fo reference.fasta
	else
		echo "Neither reference file provided nor reference genome and panel's bed file provided."
		echo "Please provide any of these files to run the pipeline."
	fi

	###################################################################
	# Edit the config file if snv_filter, reference_2 and 		  #
	# and barcode_index are provided 				  #
	###################################################################

	cat ~{config_json}

	mkdir references
	cp ~{path_to_flist} "samples.txt"
	# Add sample_id as first line of file
	sed -i '1s/^/sample_id\n/' samples.txt

	if [[ "~{reference_amplicons_2}" != '' ]]; then
		echo "Reference 2 file provided"
		python /Code/add_entry_to_json.py ~{config_json} "reference_amplicons_2" "references/reference2.fasta"
	else
		echo "Reference 2 file not provided"
	fi

	if [[ "~{path_to_snv}" != '' ]]; then
		echo "Path to SNV file provided"
		python /Code/add_entry_to_json.py ~{config_json} "path_to_snv" "references/snv_filter.tsv"
	else
		echo "Path to SNV file not provided"
	fi
	
	if [[ "~{barcodes_index}" != '' ]]; then
		echo "Barcode index file provided"
		python /Code/add_entry_to_json.py ~{config_json} "path_to_barcodes" "references/barcodes_index.csv"
	else
		echo "Barcode index file not provided"
	fi

	cat ~{config_json} >> config_json.json

	>>>

	output {
		File config_json_out = "config_json.json"
		File path_to_flist_o = "samples.txt"
		File reference_out = "reference.fasta"
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
