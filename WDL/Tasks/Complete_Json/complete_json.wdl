version 1.0

task complete_json {
	input {
		File bioinformatics_json
		Array[Boolean] use_validated_fastqs
		Array[Int] trim_galore_qvalue
		Array[Int] trim_galore_length
		Array[Int] downsample_fraction
		Int? minreads_threshold
		Int? contamination_threshold
		Int? primer_distance_threshold
		String Class
		String maxEE
		String trimRight
		Int minLen
		String truncQ
		String matchIDs
		Int max_consist
		Float omegaA
		String pool
		Int justConcatenate
		Int maxMismatch
		String strain
		String strain2
		String min_reads
		String min_samples
		String max_snv_dist
		String max_indel_dist
		Boolean include_failed
		Boolean exclude_bimeras
		String polyN

		String ampseq_export_format
		String metadata_variable1_name
		String metadata_latitude_name
		String metadata_longitude_name
		String sample_id_pat
		Int min_abd
		Float min_ratio
		String off_target_formula
		String flanking_INDEL_formula
		Int homopolymer_length
		String SNV_in_homopolymer_formula
		String INDEL_in_homopolymer_formula
		String bimera_formula
		String PCR_errors_formula
		Float sample_ampl_rate
		Float locus_ampl_rate
	}
	command <<<
	export TMPDIR=tmp
	set -euxo pipefail
	
	#####################################################
	# Make bioinformatics_info.json			    #
	#####################################################

	touch config.json
	cat "~{bioinformatics_json}" > config.json

	echo "ampseq.t_001_validate_fastqs.use_validated_fastqs" > samples_tmp.csv
	echo "~{use_validated_fastqs[1]}" >> samples_tmp.csv
	python /Code/add_entry_to_json.py samples_tmp.csv config.json

	echo "ampseq.t_002_cutadapters.trim_galore_qvalue" > samples_tmp.csv
	echo "~{trim_galore_qvalue[1]}" >> samples_tmp.csv
	python /Code/add_entry_to_json.py samples_tmp.csv config.json

	echo "ampseq.t_002_cutadapters.trim_galore_length" > samples_tmp.csv
	echo "~{trim_galore_length[1]}" >> samples_tmp.csv
	python /Code/add_entry_to_json.py samples_tmp.csv config.json

	echo "ampseq.t_002_cutadapters.trim_galore_length" > samples_tmp.csv
	echo "~{trim_galore_length[1]}" >> samples_tmp.csv
	python /Code/add_entry_to_json.py samples_tmp.csv config.json

	if grep -q "barcodes_index" config.json; then
		echo "ampseq.t_0001_contamination_detection.contamination_threshold" > samples_tmp.csv
		echo "~{contamination_threshold}" >> samples_tmp.csv
		python /Code/add_entry_to_json.py samples_tmp.csv config.json

		echo "ampseq.t_0001_contamination_detection.minreads_threshold" > samples_tmp.csv
		echo "~{minreads_threshold}" >> samples_tmp.csv
		python /Code/add_entry_to_json.py samples_tmp.csv config.json
		
		echo "ampseq.t_0001_contamination_detection.primer_distance_threshold" > samples_tmp.csv
		echo "~{primer_distance_threshold}" >> samples_tmp.csv
		python /Code/add_entry_to_json.py samples_tmp.csv config.json
	fi	

	echo "ampseq.t_004_amplicon_denoising.Class" > samples_tmp.csv
	echo "~{Class}" >> samples_tmp.csv
	python /Code/add_entry_to_json.py samples_tmp.csv config.json

	echo "ampseq.t_004_amplicon_denoising.maxEE" > samples_tmp.csv
	echo "~{maxEE}" >> samples_tmp.csv
	python /Code/add_entry_to_json.py samples_tmp.csv config.json

	echo "ampseq.t_004_amplicon_denoising.trimRight" > samples_tmp.csv
	echo "~{trimRight}" >> samples_tmp.csv
	python /Code/add_entry_to_json.py samples_tmp.csv config.json

	echo "ampseq.t_004_amplicon_denoising.minLen" > samples_tmp.csv
	echo "~{minLen}" >> samples_tmp.csv
	python /Code/add_entry_to_json.py samples_tmp.csv config.json

	echo "ampseq.t_004_amplicon_denoising.truncQ" > samples_tmp.csv
	echo "~{truncQ}" >> samples_tmp.csv
	python /Code/add_entry_to_json.py samples_tmp.csv config.json

	echo "ampseq.t_004_amplicon_denoising.pool" > samples_tmp.csv
	echo "~{pool}" >> samples_tmp.csv
	python /Code/add_entry_to_json.py samples_tmp.csv config.json

	echo "ampseq.t_004_amplicon_denoising.justConcatenate" > samples_tmp.csv
	echo "~{justConcatenate}" >> samples_tmp.csv
	python /Code/add_entry_to_json.py samples_tmp.csv config.json

	echo "ampseq.t_004_amplicon_denoising.maxMismatch" > samples_tmp.csv
	echo "~{maxMismatch}" >> samples_tmp.csv
	python /Code/add_entry_to_json.py samples_tmp.csv config.json

	echo "ampseq.t_004_amplicon_denoising.strain" > samples_tmp.csv
	echo "~{strain}" >> samples_tmp.csv
	python /Code/add_entry_to_json.py samples_tmp.csv config.json

	echo "ampseq.t_004_amplicon_denoising.strain2" > samples_tmp.csv
	echo "~{strain2}" >> samples_tmp.csv
	python /Code/add_entry_to_json.py samples_tmp.csv config.json

	echo "ampseq.t_004_amplicon_denoising.min_reads" > samples_tmp.csv
	echo "~{min_reads}" >> samples_tmp.csv
	python /Code/add_entry_to_json.py samples_tmp.csv config.json

	echo "ampseq.t_004_amplicon_denoising.min_samples" > samples_tmp.csv
	echo "~{min_samples}" >> samples_tmp.csv
	python /Code/add_entry_to_json.py samples_tmp.csv config.json

	echo "ampseq.t_004_amplicon_denoising.max_snv_dist" > samples_tmp.csv
	echo "~{max_snv_dist}" >> samples_tmp.csv
	python /Code/add_entry_to_json.py samples_tmp.csv config.json
	
	echo "ampseq.t_004_amplicon_denoising.max_indel_dist" > samples_tmp.csv
	echo "~{max_indel_dist}" >> samples_tmp.csv
	python /Code/add_entry_to_json.py samples_tmp.csv config.json

	echo "ampseq.t_004_amplicon_denoising.include_failed" > samples_tmp.csv
	echo "~{include_failed}" >> samples_tmp.csv
	python /Code/add_entry_to_json.py samples_tmp.csv config.json
	
	echo "ampseq.t_004_amplicon_denoising.exclude_bimeras" > samples_tmp.csv
	echo "~{exclude_bimeras}" >> samples_tmp.csv
	python /Code/add_entry_to_json.py samples_tmp.csv config.json

	echo "ampseq.t_004_amplicon_denoising.polyN" > samples_tmp.csv
	echo "~{polyN}" >> samples_tmp.csv
	python /Code/add_entry_to_json.py samples_tmp.csv config.json

	echo "ampseq.t_005_asv_filtering.ampseq_export_format" > samples_tmp.csv
	echo "~{ampseq_export_format}" >> samples_tmp.csv
	python /Code/add_entry_to_json.py samples_tmp.csv config.json

	echo "ampseq.t_005_asv_filtering.metadata_variable1_name" > samples_tmp.csv
	echo "~{metadata_variable1_name}" >> samples_tmp.csv
	python /Code/add_entry_to_json.py samples_tmp.csv config.json
	
	echo "ampseq.t_005_asv_filtering.metadata_latitude_name" > samples_tmp.csv
	echo "~{metadata_latitude_name}" >> samples_tmp.csv
	python /Code/add_entry_to_json.py samples_tmp.csv config.json

	echo "ampseq.t_005_asv_filtering.metadata_longitude_name" > samples_tmp.csv
	echo "~{metadata_longitude_name}" >> samples_tmp.csv
	python /Code/add_entry_to_json.py samples_tmp.csv config.json

	echo "ampseq.t_005_asv_filtering.sample_id_pat" > samples_tmp.csv
	echo "~{sample_id_pat}" >> samples_tmp.csv
	python /Code/add_entry_to_json.py samples_tmp.csv config.json

	echo "ampseq.t_005_asv_filtering.min_abd" > samples_tmp.csv
	echo "~{min_abd}" >> samples_tmp.csv
	python /Code/add_entry_to_json.py samples_tmp.csv config.json

	echo "ampseq.t_005_asv_filtering.min_ratio" > samples_tmp.csv
	echo "~{min_ratio}" >> samples_tmp.csv
	python /Code/add_entry_to_json.py samples_tmp.csv config.json

	echo "ampseq.t_005_asv_filtering.off_target_formula" > samples_tmp.csv
	echo "~{off_target_formula}" >> samples_tmp.csv
	python /Code/add_entry_to_json.py samples_tmp.csv config.json

	echo "ampseq.t_005_asv_filtering.flanking_INDEL_formula" > samples_tmp.csv
	echo "~{flanking_INDEL_formula}" >> samples_tmp.csv
	python /Code/add_entry_to_json.py samples_tmp.csv config.json

	echo "ampseq.t_005_asv_filtering.homopolymer_length" > samples_tmp.csv
	echo "~{homopolymer_length}" >> samples_tmp.csv
	python /Code/add_entry_to_json.py samples_tmp.csv config.json

	echo "ampseq.t_005_asv_filtering.SNV_in_homopolymer_formula" > samples_tmp.csv
	echo "~{SNV_in_homopolymer_formula}" >> samples_tmp.csv
	python /Code/add_entry_to_json.py samples_tmp.csv config.json

	echo "ampseq.t_005_asv_filtering.INDEL_in_homopolymer_formula" > samples_tmp.csv
	echo "~{INDEL_in_homopolymer_formula}" >> samples_tmp.csv
	python /Code/add_entry_to_json.py samples_tmp.csv config.json

	echo "ampseq.t_005_asv_filtering.bimera_formula" > samples_tmp.csv
	echo "~{bimera_formula}" >> samples_tmp.csv
	python /Code/add_entry_to_json.py samples_tmp.csv config.json

	echo "ampseq.t_005_asv_filtering.PCR_errors_formula" > samples_tmp.csv
	echo "~{PCR_errors_formula}" >> samples_tmp.csv
	python /Code/add_entry_to_json.py samples_tmp.csv config.json

	echo "ampseq.t_005_asv_filtering.sample_ampl_rate" > samples_tmp.csv
	echo "~{sample_ampl_rate}" >> samples_tmp.csv
	python /Code/add_entry_to_json.py samples_tmp.csv config.json

	echo "ampseq.t_005_asv_filtering.locus_ampl_rate" > samples_tmp.csv
	echo "~{locus_ampl_rate}" >> samples_tmp.csv
	python /Code/add_entry_to_json.py samples_tmp.csv config.json

	mv config.json bioinformatics_info.json

	>>>

	output {
		File bioinformatics_info_json_o = "bioinformatics_info.json"
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

