version 1.0

workflow ampseq {
	input {	
		#General commands
		Array[File] path_to_r1
		Array[File] path_to_r2
		File path_to_flist
		File pr1
		File pr2
		File reference1
		File reference2
		File path_to_snv
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

		#Variables to reconstruct data table
		Array[String] Sample_id
		Array[String] Geo_Level
		Array[String] Temp_Level
		Array[String] Longitude
		Array[String] Latitude
	}

	call ampseq_pipeline {
		input:
			path_to_r1 = path_to_r1,
			path_to_r2 = path_to_r2,
			path_to_flist = path_to_flist,
			pr1 = pr1,
			pr2 = pr2,
			reference1 = reference1,
			reference2 = reference2,
			path_to_snv = path_to_snv,
			run_id = run_id,
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
	
	call make_metadata {
		input:			
			run_id = run_id,
			Sample_id = Sample_id,
			Geo_Level = Geo_Level,
			Temp_Level = Temp_Level,
			Longitude = Longitude,
			Latitude = Latitude
	}

	output {
		File ASVBimeras_f = ampseq_pipeline.ASVBimeras
		File CIGARVariants_Bfilter_f = ampseq_pipeline.CIGARVariants_Bfilter
		File ASV_to_CIGAR_f = ampseq_pipeline.ASV_to_CIGAR
		File seqtab_f = ampseq_pipeline.seqtab
		File ASVTable_f = ampseq_pipeline.ASVTable
		File ASVSeqs_f = ampseq_pipeline.ASVSeqs
		File missing_files_f = ampseq_pipeline.missing_files
		File? decontamination_sample_cards_f = ampseq_pipeline.decontamination_sample_cards
		File? decontamination_report_f = ampseq_pipeline.decontamination_report
		File? metadata_f = make_metadata.metadata_out
	}
}

task ampseq_pipeline {
	input {
		Array[File] path_to_r1
		Array[File] path_to_r2
		File path_to_flist
		File pr1
		File pr2
		File reference1
		File reference2
		File path_to_snv
		Array[String] run_id
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
		"path_to_flist": sub(path_to_flist, "gs://", "/cromwell_root/"),
		"pr1": sub(pr1, "gs://", "/cromwell_root/"),
		"pr2": sub(pr2, "gs://", "/cromwell_root/"),
		"reference1": sub(reference1, "gs://", "/cromwell_root/"),
		"reference2": sub(reference2, "gs://", "/cromwell_root/"),
		"path_to_snv": sub(path_to_snv, "gs://", "/cromwell_root/"),
#		"path_to_flist": "barcodes_matches.csv", 
#		"pr1": "primers_fw.fasta", 
#		"pr2": "primers_rv.fasta", 
#		"reference1": "reference1.fasta", 
#		"reference2": "reference2.fasta", 
#		"path_to_snv": "path_to_snv.fasta", 
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
	mkdir fq_dir

	gsutil -m cp -r ~{sep = ' ' path_to_r1} fq_dir/
	gsutil -m cp -r ~{sep = ' ' path_to_r2} fq_dir/

	#cp ~{sep = ' ' path_to_r1} fq_dir/
	#cp ~{sep = ' ' path_to_r2} fq_dir/
	#cat ~{config_json}
	#echo $PWD
	#cp ~{path_to_flist} barcodes_matches.csv
	#cp ~{pr1} primers_fw.fasta
	#cp ~{pr2} primers_rv.fasta
	#cp ~{reference1} reference1.fasta
	#cp ~{reference2} reference2.fasta
	#cp ~{path_to_snv} path_to_snv.fasta
	#ls 

	#Move reference files to the main level	
	#Check if the first line in barcodes_matches.csv indicates the presence of inline barcodes
	if grep -q "," ~{path_to_flist} ; then
		echo "Sequencing run with inline barcodes. Performing analysis of combinatorial indices followed by denoising"
		find . -type f
		python /Code/Amplicon_TerraPipeline.py --config ~{config_json} --terra --meta --adaptor_removal --contamination --separate_reads --primer_removal --dada2 --postproc_dada2 --asv_to_cigar
		find . -type f
		Rscript /Code/render_report.R -d /cromwell_root/Report/Merge/ -o /cromwell_root/Report/ -p ~{path_to_flist} -m 1000 -c 0.5 -mf /cromwell_root/Results/missing_files.tsv
		find . -type f
		tar -czvf Report_Cards.tar.gz /cromwell_root/Report
		find . -type f
	else
		echo "Sequencing run without inline barcodes. Skipping analysis of combinatorial indices and performing only denoising"
		find . -type f
		python /Code/Amplicon_TerraPipeline.py --config ~{config_json} --terra --meta --adaptor_removal --primer_removal --dada2 --postproc_dada2 --asv_to_cigar
		#python /Code/Amplicon_TerraPipeline.py --config ~{config_json} --terra --meta --adaptor_removal --separate_reads --primer_removal --dada2 --postproc_dada2 --asv_to_cigar
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
		File missing_files = "Results/missing_files.tsv"
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

task make_metadata {
	input {
		#Variables to reconstruct data table
		Array[String] Sample_id
		Array[String] Geo_Level
		Array[String] Temp_Level
		Array[String] Longitude
		Array[String] Latitude
		Array[String] run_id
	}

	Array[Array[String]] metadata_array = transpose([Sample_id, Geo_Level, Temp_Level, Longitude, Latitude])
        File metadata = write_tsv(metadata_array)

	command <<<
		sed 's/\t/,/g' ~{metadata} > metadata.csv

		run_id_array=(~{sep = ' ' run_id})
		unique_id=$(printf "%s\n" "${run_id_array[@]}" | sort -u | tr '\n' '_')
		unique_id="${unique_id%_}"
	
		cp metadata.csv "${unique_id}_metadata.csv"
	>>>
	output {
 		File metadata_out = glob("*.csv")[0]
	}
	runtime {
		docker: "broadinstitute/horsefish"
	}
}
