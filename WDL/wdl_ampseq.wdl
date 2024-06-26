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

	if(run_demultiplexing) {
		call ampseq_pipeline_demult {
			input: 
				config_json = prepare_files.config_json_out,
				path_to_r1 = path_to_r1,
				path_to_r2 = path_to_r2,
				pr1 = pr1,
				pr2 = pr2,
				reference = prepare_files.reference_out,
				reference2 = reference2,
				path_to_snv = path_to_snv
		}
	}

	if(!run_demultiplexing) {
		call ampseq_pipeline_no_demult {
			input:
				config_json = prepare_files.config_json_out,
				path_to_r1 = path_to_r1,
				path_to_r2 = path_to_r2,
				pr1 = pr1,
				pr2 = pr2,
				reference = prepare_files.reference_out,
				reference2 = reference2,
				path_to_snv = path_to_snv
		}
	}

	call ampseq_pipeline_denoise {
		input:
			config_json = prepare_files.config_json_out,
			path_to_r1 = path_to_r1,
			path_to_r2 = path_to_r2,
			run_id = run_id,
			primer_rem = if (run_demultiplexing) then ampseq_pipeline_demult.PrimerRem else ampseq_pipeline_no_demult.PrimerRem
	}
	
	output {
		File ASVBimeras_f = ampseq_pipeline_denoise.ASVBimeras
		File CIGARVariants_Bfilter_f = ampseq_pipeline_denoise.CIGARVariants_Bfilter
		File ASV_to_CIGAR_f = ampseq_pipeline_denoise.ASV_to_CIGAR
		File seqtab_f = ampseq_pipeline_denoise.seqtab
		File ASVTable_f = ampseq_pipeline_denoise.ASVTable
		File ASVSeqs_f = ampseq_pipeline_denoise.ASVSeqs
		File? missing_files_f = if (run_demultiplexing) then ampseq_pipeline_demult.missing_files else ampseq_pipeline_no_demult.missing_files
		File? decontamination_sample_cards_f = if (run_demultiplexing) then ampseq_pipeline_demult.decontamination_sample_cards else ampseq_pipeline_no_demult.decontamination_sample_cards
		File? decontamination_report_f = if (run_demultiplexing) then ampseq_pipeline_demult.decontamination_report else ampseq_pipeline_no_demult.decontamination_report
		File? panel_reference_fasta_f = prepare_files.reference_out
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

	Map[String, String] in_map = {
		"path_to_fq": "fq_dir",
		"path_to_flist": "references/barcodes_matches.csv",
		"pr1": "references/pr1.fasta",
		"pr2": "references/pr2.fasta",
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
	##Make reference fasta file if reference not provided by user    ##
	###################################################################

	if [[ "~{reference}" != '' ]]; then
		echo "Reference file not provided."
		echo "Reference genome and panel's bed file provided."
		echo "Creating reference file"
		bedtools getfasta -fi ~{reference_genome} -bed ~{panel_bedfile} -fo reference.fasta
	elif [[ "~{reference}" == '' && "~{reference_genome}" != '' && "~{panel_bedfile}" != '' ]]; then
		echo "Reference file provided."
		cp ~{reference} reference.fasta
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

	cat ~{config_json}

	>>>

	output {
		File config_json_out = config_json
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

	./Code/file_staging.sh

	echo "Demultiplexing not requested."
	echo "No demultiplexing will be performed in the data. Read pairs assumed to be long enough to overlap and be merged."

	#Check if the first line in barcodes_matches.csv indicates the presence of inline barcodes
	if grep -q "," references/barcodes_matches ; then
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

task ampseq_pipeline_demult {
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

	./Code/file_staging.sh

	echo "Demultiplexing requested."
	echo "Demultiplexing will be performed in the data. Some read pairs assumed to be too short to overlap and be merged."

	#Check if the first line in barcodes_matches.csv indicates the presence of inline barcodes
	if grep -q "," references/barcodes_matches ; then
		echo "Sequencing run with inline barcodes. Performing analysis of combinatorial indices followed by denoising"
		find . -type f
		python /Code/Amplicon_TerraPipeline.py --config ~{config_json} --terra --meta --adaptor_removal --contamination --demultiplex --primer_removal
		find . -type f
		Rscript /Code/render_report.R -d /cromwell_root/Report/Merge/ -o /cromwell_root/Report/ -p ~{path_to_flist} -m 1000 -c 0.5 -mf /cromwell_root/Results/missing_files.tsv
		find . -type f
		tar -czvf Report_Cards.tar.gz /cromwell_root/Report
		find . -type f
	else
		echo "Sequencing run without inline barcodes. Skipping analysis of combinatorial indices and performing only denoising"
		find . -type f
		python /Code/Amplicon_TerraPipeline.py --config ~{config_json} --terra --meta --adaptor_removal --demultiplex --primer_removal
		find . -type f
	fi
	
	>>>
	output {
		File missing_files = "Results/missing_files.tsv"
		Array[File] PrimerRem = glob("Results/PrimerRem/*")
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

task ampseq_pipeline_denoise {
	input {
		Array[File] path_to_r1
		Array[File] path_to_r2
		Array[File]? primer_rem
		Array[String] run_id
		File config_json
	}

	command <<<
	set -euxo pipefail
	cat ~{config_json}

	###################################################################
	##Copy files to the working directory and run the AmpSeq pipeline##
	###################################################################

	mkdir -p Results/PrimerRem
	mkdir fq_dir
	gsutil -m cp -r ~{sep = ' ' path_to_r1} fq_dir/
	gsutil -m cp -r ~{sep = ' ' path_to_r2} fq_dir/
	gsutil -m cp -r ~{sep = ' ' primer_rem} Results/PrimerRem/

	if [ -e "Results/PrimerRem/mixed_nop_prim_meta.tsv" ];then
		echo "Demultiplex performed in the data. Some read pairs assumed to be too short to overlap and be merged."
		find . -type f
		python /Code/Amplicon_TerraPipeline.py --config ~{config_json} --terra --dada2 --demultiplexed --postproc_dada2 --asv_to_cigar
		find . -type f
	else
		echo "Demultiplex not performed in the data. Read pairs assumed to be long enough to overlap and be merged."
		find . -type f
		python /Code/Amplicon_TerraPipeline.py --config ~{config_json} --terra --dada2 --postproc_dada2 --asv_to_cigar
		find . -type f
	fi

	run_id_array=(~{sep = ' ' run_id})
	unique_id=$(printf "%s\n" "${run_id_array[@]}" | sort -u | tr '\n' '_')
	unique_id="${unique_id%_}"
	
	cp Results/CIGARVariants_Bfilter.out.tsv "${unique_id}_CIGARVariants_Bfilter.out.tsv"

	>>>

	output {
		File ASVBimeras = "Results/ASVBimeras.txt"
		File CIGARVariants_Bfilter = glob("*.out.tsv")[0]
		File ASV_to_CIGAR = "Results/ASV_to_CIGAR/ASV_to_CIGAR.out.txt"
		File seqtab = "Results/seqtab.tsv"
		File ASVTable = "Results/PostProc_DADA2/ASVTable.txt"
		File ASVSeqs = "Results/PostProc_DADA2/ASVSeqs.fasta"
		File missing_files = missing_files
		File? decontamination_sample_cards = decontamination_sample_cards
		File? decontamination_report = decontamination_report
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

