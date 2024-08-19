version 1.0

workflow ampseq {
	input {	
		#General commands
		Array[File] path_to_r1
		Array[File] path_to_r2
		File path_to_flist
		File pr1
		File pr2
		File? reference
		File? reference2
		File? path_to_snv
		File? reference_genome
		File? panel_bedfile
		Array[String] run_id

		String pattern_fw = "*_L001_R1_001.fastq.gz"
		String pattern_rv = "*_L001_R2_001.fastq.gz"

		#Commands for AmpSeq
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
		String adapter = "None"

		#Command for the decontamination pipeline
		Int minreads_threshold = 1000
		Float contamination_threshold = 0.5
		String verbose = "False"
		Boolean run_demultiplexing = false
	}

	call prepare_files {
		input:
			path_to_r1 = path_to_r1,
			path_to_r2 = path_to_r2,
			path_to_flist = path_to_flist,
			pr1 = pr1,
			pr2 = pr2,
			reference = reference,
			reference2 = reference2,
			path_to_snv = path_to_snv,
			reference_genome = reference_genome,
			panel_bedfile = panel_bedfile,
			pattern_fw = pattern_fw,
			pattern_rv = pattern_rv,
			Class = Class,
			maxEE = maxEE,
			trimRight = trimRight,
			minLen = minLen,
			truncQ = truncQ,
			matchIDs = matchIDs,
			max_consist = max_consist,
			omegaA = omegaA,
			saveRdata = saveRdata,
			justConcatenate = justConcatenate,
			maxMismatch = maxMismatch,
			no_ref = no_ref,
			adjust_mode = adjust_mode,
			strain = strain,
			strain2 = strain2,
			polyN = polyN,
			min_reads = min_reads,
			min_samples = min_samples,
			max_snv_dist = max_snv_dist,
			max_indel_dist = max_indel_dist,
			include_failed = include_failed,
			exclude_bimeras = exclude_bimeras,
			minreads_threshold = minreads_threshold,
			contamination_threshold = contamination_threshold,
			verbose = verbose,
			adapter = adapter
	}

	if(!run_demultiplexing) {
		call ampseq_pipeline_no_demult {
			input:
				config_json = prepare_files.config_json_out,
				path_to_flist = path_to_flist,
				path_to_r1 = path_to_r1,
				path_to_r2 = path_to_r2,
				pr1 = pr1,
				pr2 = pr2,
				reference = prepare_files.reference_out,
				reference2 = reference2,
				path_to_snv = path_to_snv
		}
	}

#	if(run_demultiplexing) {
#		call ampseq_pipeline_demult {
#			input: 
#				config_json = prepare_files.config_json_out,
#				path_to_flist = path_to_flist,
#				path_to_r1 = path_to_r1,
#				path_to_r2 = path_to_r2,
#				pr1 = pr1,
#				pr2 = pr2,
#				reference = prepare_files.reference_out,
#				reference2 = reference2,
#				path_to_snv = path_to_snv
#		}
#	}
##
	call ampseq_pipeline_denoise {
		input:
			config_json = prepare_files.config_json_out,
			path_to_flist = path_to_flist,
			path_to_r1 = path_to_r1,
			path_to_r2 = path_to_r2,
			pr1 = pr1,
			pr2 = pr2,
			reference = prepare_files.reference_out,
			reference2 = reference2,
			run_id = run_id,
			path_to_snv = path_to_snv,
	#		primer_rem = if (run_demultiplexing) then ampseq_pipeline_demult.PrimerRem else ampseq_pipeline_no_demult.PrimerRem,
	#		adaptor_rem = if (run_demultiplexing) then ampseq_pipeline_demult.AdaptorRem else ampseq_pipeline_no_demult.AdaptorRem

			###REMOVE THIS VARIABLES AFTER TESTING###

			primer_rem = ampseq_pipeline_no_demult.PrimerRem,
			adaptor_rem = ampseq_pipeline_no_demult.AdaptorRem
	}

	output {
		File? panel_reference_fasta_f = prepare_files.reference_out
		
		# DADA2
		File ASVBimeras_f = ampseq_pipeline_denoise.ASVBimeras
		File seqtab_f = ampseq_pipeline_denoise.seqtab

		# PostProc_DADA2
		File ASVTable_f = ampseq_pipeline_denoise.ASVTable
		File ASVSeqs_f = ampseq_pipeline_denoise.ASVSeqs

		# ASV_to_CIGAR
		File CIGARVariants_Bfilter_f = ampseq_pipeline_denoise.CIGARVariants_Bfilter
		File ASV_to_CIGAR_f = ampseq_pipeline_denoise.ASV_to_CIGAR
		File ZeroReadsSampleList_f = ampseq_pipeline_denoise.ZeroReadsSampleList

		###REMOVE THIS VARIABLES AFTER TESTING###
		File config_json_out_f = prepare_files.config_json_out
		File? missing_files_f = ampseq_pipeline_no_demult.missing_files
		File? decontamination_sample_cards_f = ampseq_pipeline_no_demult.decontamination_sample_cards
		File? decontamination_report_f = ampseq_pipeline_no_demult.decontamination_report
	}
}

task prepare_files {
	input {
		Array[File] path_to_r1
		Array[File] path_to_r2
		File path_to_flist
		File pr1
		File pr2
		File? reference
		File? reference2
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
	}

	String flist_refdir = "references/" + basename(path_to_flist)
	String pr1_refdir = "references/" + basename(pr1)
	String pr2_refdir = "references/" + basename(pr2)

	Map[String, String] in_map = {
		"path_to_fq": "fq_dir",
		"path_to_flist": flist_refdir,
		"pr1": pr1_refdir,
		"pr2": pr2_refdir,
		"reference": "references/reference.fasta",
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
	if has_bom ~{panel_bedfile}; then 
		echo "File has a BOM. Removing."
		sed -i 's/\xef\xbb\xbf//' ~{panel_bedfile}
	else
		echo "File does not have a BOM. Skipping this step."
	fi

	###################################################################
	##Make reference fasta file if reference not provided by user    ##
	###################################################################

	if [[ "~{reference}" != '' ]]; then
		echo "Reference file provided."
		cp ~{reference} reference.fasta
	elif [[ "~{reference}" == '' && "~{reference_genome}" != '' && "~{panel_bedfile}" != '' ]]; then
		echo "Reference file not provided."
		echo "Reference genome and panel's bed file provided."
		echo "Creating reference file"
		bedtools getfasta -fi ~{reference_genome} -bed ~{panel_bedfile} -fo reference.fasta
	else
		echo "Neither reference file provided nor reference genome and panel's bed file provided."
		echo "Please provide any of these files to run the pipeline."
	fi

	###################################################################
	##Edit the config file if snv_filter and reference_2 are provided##
	###################################################################

	cat ~{config_json}

	if [[ "~{reference2}" != '' ]]; then
		echo "Reference 2 file provided"
		python /Code/add_entry_to_json.py ~{config_json} "reference2" "references/reference2.fasta"
	else
		echo "Reference 2 file not provided"
	fi

	if [[ "~{path_to_snv}" != '' ]]; then
		echo "Path to SNV file provided"
		python /Code/add_entry_to_json.py ~{config_json} "path_to_snv" "references/snv_filter.tsv"
	else
		echo "Path to SNV file not provided"
	fi

	cat ~{config_json} >> config_json.json

	>>>

	output {
		File config_json_out = "config_json.json"
		File reference_out = "reference.fasta"
	}

	runtime {
		cpu: 1
		memory: "1 GB"
		disks: "local-disk 10 HDD"
		bootDiskSizeGb: 10
		preemptible: 3
		maxRetries: 1
		docker: "jorgeamaya/fileprep_ampseq"
	}
}

task ampseq_pipeline_no_demult {
	input {
		Array[File] path_to_r1
		Array[File] path_to_r2
		File path_to_flist
		File pr1
		File pr2
		File reference
		File? reference2
		File? path_to_snv
		File config_json
	}

	command <<<
	set -euxo pipefail
	cat ~{config_json}

	mkdir fq_dir
	mkdir references
	gsutil -m cp -r ~{sep = ' ' path_to_r1} fq_dir/
	gsutil -m cp -r ~{sep = ' ' path_to_r2} fq_dir/
	gsutil cp ~{path_to_flist} references/
	gsutil cp ~{pr1} references/
	gsutil -m cp -r ~{pr2} references/
	gsutil -m cp -r ~{reference} references/

	~{"gsutil cp " + reference2 + " references/"}
	~{"gsutil cp " + path_to_snv + " references/"}

	echo "Demultiplexing not requested."
	echo "No demultiplexing will be performed in the data. Read pairs assumed to be long enough to overlap and be merged."

	#Check if the first line in barcodes_matches.csv indicates the presence of inline barcodes
	if grep -q "," "references/$(basename -- ~{path_to_flist})" ; then
		echo "Sequencing run with inline barcodes. Performing analysis of combinatorial indices."
		find . -type f
		python /Code/Amplicon_TerraPipeline.py --config ~{config_json} --terra --meta --adaptor_removal --contamination --primer_removal
		find . -type f
		Rscript /Code/render_report.R -d /cromwell_root/Report/Merge/ -o /cromwell_root/Report/ -p ~{path_to_flist} -m 1000 -c 0.5 -mf /cromwell_root/Results/missing_files.tsv
		find . -type f
		tar -czvf Report_Cards.tar.gz /cromwell_root/Report
		find . -type f
	else
		echo "Sequencing run without inline barcodes. Skipping analysis of combinatorial indices."
		find . -type f
		python /Code/Amplicon_TerraPipeline.py --config ~{config_json} --terra --meta --adaptor_removal --primer_removal
		find . -type f
	fi
	
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
		preemptible: 3
		maxRetries: 1
		docker: 'jorgeamaya/ampseq'
	}
}

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
#		docker: 'jorgeamaya/ampseq'
#	}
#}
#
task ampseq_pipeline_denoise {
	input {
		Array[File] path_to_r1
		Array[File] path_to_r2
		File path_to_flist
		File pr1
		File pr2
		File reference
		File? reference2
		File? path_to_snv
		Array[File]? primer_rem
		Array[File]? adaptor_rem
		Array[String] run_id
		File config_json
	}

	command <<<
	set -euxo pipefail
	cat ~{config_json}

	###################################################################
	##Copy files to the working directory and run the AmpSeq pipeline##
	###################################################################
	mkdir -p Results/
	mkdir -p Results/PrimerRem
	mkdir -p Results/AdaptorRem
	mkdir -p Results/PostProc_DADA2
	mkdir -p Results/ASV_to_CIGAR
	mkdir -p Results/alignments

	mkdir fq_dir
	mkdir references
	gsutil -m cp -r ~{sep = ' ' path_to_r1} fq_dir/
	gsutil -m cp -r ~{sep = ' ' path_to_r2} fq_dir/
	gsutil cp ~{path_to_flist} references/
	gsutil cp ~{pr1} references/
	gsutil cp ~{pr2} references/
	gsutil -m cp -r ~{sep = ' ' primer_rem} Results/PrimerRem/
	gsutil -m cp -r ~{sep = ' ' adaptor_rem} Results/AdaptorRem/
	gsutil cp ~{reference} references/

	~{"gsutil cp " + reference2 + " references/"}
	~{"gsutil cp " + path_to_snv + " references/"}

	if [ -e "Results/PrimerRem/mixed_nop_prim_meta.tsv" ]; then
		echo "Demultiplex performed in the data. Some read pairs assumed to be too short to overlap and be merged."
		python /Code/Amplicon_TerraPipeline.py --config ~{config_json} --terra --dada2 --demultiplexed
	else
		echo "Demultiplex not performed in the data. Read pairs assumed to be long enough to overlap and be merged."
		python /Code/Amplicon_TerraPipeline.py --config ~{config_json} --terra --dada2 
	fi

	# Run postproc_DADA2
	echo "Running post processing of DADA2 results..."
	python /Code/Amplicon_TerraPipeline.py --config ~{config_json} --terra --postproc_dada2

	# Run ASV_to_CIGAR
	echo "Converting ASV to CIGAR tables..."
	python /Code/Amplicon_TerraPipeline.py --config ~{config_json} --terra --asv_to_cigar

	#run_id_array=(~{sep = ' ' run_id})
	#unique_id=$(printf "%s\n" "${run_id_array[@]}" | sort -u | tr '\n' '_')
	#unique_id="${unique_id%_}"
	#cp Results/CIGARVariants_Bfilter.out.tsv "${unique_id}_CIGARVariants_Bfilter.out.tsv"

	>>>

	output {
		File ASVBimeras = "Results/ASVBimeras.txt"
		File CIGARVariants_Bfilter = glob("*.out.tsv")[0]
		File ASV_to_CIGAR = "Results/ASV_to_CIGAR/ASV_to_CIGAR.out.txt"
		File ZeroReadsSampleList = "Results/ASV_to_CIGAR/ZeroReadsSampleList.txt"

		File seqtab = "Results/seqtab.tsv"
		File ASVTable = "Results/PostProc_DADA2/ASVTable.txt"
		File ASVSeqs = "Results/PostProc_DADA2/ASVSeqs.fasta"
	}

	runtime {
		cpu: 1
		memory: "15 GiB"
		disks: "local-disk 10 HDD"
		bootDiskSizeGb: 10
		preemptible: 3
		maxRetries: 1
		docker: 'jorgeamaya/ampseq'
	}
}

task ampseq_pipeline_asv_filtering {
	input {
		String out_prefix 
		File? panel_bedfile
		File reference		#[TODO: Ask about compatibility for second reference panel (i.e. reference2)]

		# Results of post-processing and CIGAR conversion
		File CIGARVariants
		File ASVTable
		File ASV_to_CIGAR
		File ASVSeqs
		File ZeroReadsSampleList

		# MHap ASV filtering thresholds [TODO: Ask if defaults are appropriate]
		String sample_id_pat = '.'
		Int? min_abd = 10
		Float? min_ratio = 0.1
		String? off_target_formula = "dVSITES_ij>=0.3"
		String? flanking_INDEL_formula = "flanking_INDEL==TRUE&h_ij>=0.66"
		Int? homopolymer_length = 5
		String? SNV_in_homopolymer_formula = "SNV_in_homopolymer==TRUE&h_ij>=0.66"
		String? INDEL_in_homopolymer_formula = "INDEL_in_homopolymer==TRUE&h_ij>=0.66"
		String? bimera_formula = "bimera==TRUE&h_ij>=0.66"
		String? PCR_errors_formula = "h_ij>=0.66&h_ijminor>=0.66&p_ij>=0.05"
		Float? sample_ampl_rate = 0.75
		Float? locus_ampl_rate = 0.75
	}

	File ref_for_markers = select_first([panel_bedfile, reference])
	###########################################
	# MHap - Define appropriate directories
	String wd = "Results/"
	String fd = "Code/MHap/"
	String rd = "references/"
	# The directories below are subdirectories of "wd" (from MHap specs)
	String cigar_variants_dir = "cigar_variants/"
	String asv_table_dir = "asv_tables/"
	String asv2cigar_dir = "asv2cigar/"
	String asv_seq_dir = "asv_seq/"
	String zero_read_sample_list_dir = "zeroReadSampleList/"
	# Output filepath
	File ampseq_obj_path = "Results/" + out_prefix + ".csv"
	###########################################

	# TO-DO: Check if selecting such defaults are appropriate.

	command <<<
		set -euxo pipefail

		# Create directories - required for MHap script
		mkdir -p Results/
		mkdir -p Results/~{cigar_variants_dir}
		mkdir -p Results/~{asv_table_dir}
		mkdir -p Results/~{asv2cigar_dir}
		mkdir -p Results/~{asv_seq_dir}
		mkdir -p Results/~{zero_read_sample_list_dir}
		mkdir -p references/

		# Copy correct files to proper directories - required for MHap script
		gsutil cp ~{reference} references/
		~{"gsutil cp " + panel_bedfile + " references/"}

		gsutil cp ~{CIGARVariants} Results/~{cigar_variants_dir}
		gsutil cp ~{ASVTable} Results/~{ASVTable}
		gsutil cp ~{ASV_to_CIGAR} Results/~{asv2cigar_dir}
		gsutil cp ~{ASVSeqs} Results/~{asv_seq_dir}
		gsutil cp ~{ZeroReadsSampleList} Results/~{zero_read_sample_list_dir}

		# Create marker table
		python /Code/markersTable_from_bed.py -i ~{ref_for_markers} -o Results/markersTable.csv

		#Call MHap_Analysis_pipeline
		echo "Applying filters to ASVs..."
		Rscript /Code/MHap_Analysis_pipeline.R \
		-fd ~{fd} \
		-wd ~{wd} \
		-rd ~{rd} \
		-cigar_files ~{cigar_variants_dir} \
		-asv_table_files ~{asv_table_dir} \
		-asv2cigar_files ~{asv2cigar_dir} \
		-asv_seq_files ~{asv_seq_dir} \
		-zero_read_sample_list ~{zero_read_sample_list_dir} \
		-o ~{out_prefix} \
		-markers Results/markersTable.csv \
		-sample_id_pattern ~{sample_id_pat} \
		~{"-min_abd " + min_abd} \
		~{"-min_ratio " + min_ratio} \
		~{"-off_target_formula " + off_target_formula} \
		~{"-flanking_INDEL_formula " + flanking_INDEL_formula} \
		~{"-homopolymer_length " + homopolymer_length} \
		~{"-SNV_in_homopolymer_formula " + SNV_in_homopolymer_formula} \
		~{"-INDEL_in_homopolymer_formula " + INDEL_in_homopolymer_formula} \
		~{"-bimera_formula " + bimera_formula} \
		~{"-PCR_errors_formula " + PCR_errors_formula} \
		~{"-samprate " + sample_ampl_rate} \
		~{"-lamprate " + locus_ampl_rate}

		echo "Finished filtering ASVs!"
		
	>>>
	output {
		File markersTable = "Results/markersTable.csv"
		File ampseq_object = ampseq_obj_path
	}
	runtime {
		cpu: 1
		memory: "40 GiB"
		disks: "local-disk 10 HDD"
		bootDiskSizeGb: 10
		preemptible: 3
		maxRetries: 1
		docker: 'jorgeamaya/ampseq'
	}
}