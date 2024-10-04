version 1.0

task asv_filtering {
	input {
		String out_prefix 
		File? panel_bedfile
		File? reference		#[TODO: Ask about compatibility for second reference panel (i.e. reference2)]
		File? markersTable
		File reference_genome
		String? ampseq_export_format = 'excel'

		# Metadata columns	
		File sample_metadata
		String? metadata_variable1_name = 'Country'
		String? metadata_variable2_name
		String? metadata_latitude_name = 'Latitude'
		String? metadata_longitude_name = 'Longitude'

		# Results of post-processing and CIGAR conversion
		File CIGARVariants
		File ASVTable
		File ASV_to_CIGAR
		File ASVSeqs
		File ZeroReadsSampleList

		# MHap ASV filtering thresholds 
		String sample_id_pat = '.'
		Int? min_abd = 10
		Float? min_ratio = 0.1
		String? off_target_formula = "dVSITES_ij>=0.3"
		String? flanking_INDEL_formula = "flanking_INDEL==TRUE\&h_ij>=0.66"
		Int? homopolymer_length = 5
		String? SNV_in_homopolymer_formula = "SNV_in_homopolymer==TRUE\&h_ij>=0.66"
		String? INDEL_in_homopolymer_formula = "INDEL_in_homopolymer==TRUE\&h_ij>=0.66"
		String? bimera_formula = "bimera==TRUE&h_ij>=0.66"
		String? PCR_errors_formula = "h_ij>=0.66\&h_ijminor>=0.66\&p_ij>=0.05"

		#[TODO: Migrate these filters to MHap pipeline]
		Float? sample_ampl_rate = 0.3
		Float? locus_ampl_rate = 0.3
	}

	File ref_for_markers = select_first([panel_bedfile, reference])
	###########################################
	# MHap - Define appropriate directories
	String wd = "Results/"
	String fd = "/Code/MHap"
	String rd = "references/"
	String ref_genome_base = basename(reference_genome)
	String ref_base = if defined(reference) then basename(select_first([reference])) else ""
	# The directories below are subdirectories of "wd" (from MHap specs)
	String cigar_variants_dir = "cigar_variants/"
	String asv_table_dir = "asv_tables/"
	String asv2cigar_dir = "asv2cigar/"
	String asv_seq_dir = "asv_seq/"
	String zero_read_sample_list_dir = "zeroReadSampleList/"

	###########################################
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
		gsutil cp ~{reference_genome} references/
		~{"gsutil cp " + panel_bedfile + " references/"}
		gsutil cp ~{sample_metadata} Results/metadata.csv
		

		gsutil cp ~{CIGARVariants} Results/~{cigar_variants_dir}/~{out_prefix}_CIGARVariants_Bfilter.out.tsv
		gsutil cp ~{ASVTable} Results/~{asv_table_dir}
		gsutil cp ~{ASV_to_CIGAR} Results/~{asv2cigar_dir}
		gsutil cp ~{ASVSeqs} Results/~{asv_seq_dir}
		gsutil cp ~{ZeroReadsSampleList} Results/~{zero_read_sample_list_dir}
		~{"gsutil cp " + markersTable + " references/markersTable.csv"}

		# Create marker table
		if [ -f "references/markersTable.csv" ]; then
			echo "Markers table provided! Skipping creation of marker table from reference..."
		else
			echo "Markers table not provided. Creating marker table from reference..."
			python /Code/createMarkersTable.py -i ~{ref_for_markers} -o references/markersTable.csv
			echo "Finished creating markers table from reference."
		fi

		# Call MHap_Analysis_pipeline from Amplicon_TerraPipeline
		find . -type f
		find /cromwell_root -type f
		ls -l1 /cromwell_root
		root_dir=$(pwd)
		echo "Applying filters to ASVs..."
		Rscript /Code/MHap_Analysis_pipeline.R \
		-fd "~{fd}" \
		-wd "$root_dir/~{wd}" \
		-rd "$root_dir/~{rd}" \
		-cigar_files '~{cigar_variants_dir}' \
		-asv_table_files '~{asv_table_dir}' \
		-asv2cigar_files '~{asv2cigar_dir}' \
		-asv_seq_files '~{asv_seq_dir}' \
		-zero_read_sample_list '~{zero_read_sample_list_dir}' \
		-o '~{out_prefix}' \
		-markers markersTable.csv \
		-sample_id_pattern '~{sample_id_pat}' \
		~{"-min_abd " + min_abd} \
		~{"-min_ratio " + min_ratio} \
		~{"-off_target_formula \'" + off_target_formula + "\'"} \
		~{"-flanking_INDEL_formula \'" + flanking_INDEL_formula + "\'"} \
		~{"-homopolymer_length \'" + homopolymer_length + "\'"} \
		~{"-SNV_in_homopolymer_formula \'" + SNV_in_homopolymer_formula + "\'"} \
		~{"-INDEL_in_homopolymer_formula \'" + INDEL_in_homopolymer_formula + "\'"} \
		~{"-bimera_formula \'" + bimera_formula + "\'"} \
		~{"-PCR_errors_formula \'" + PCR_errors_formula + "\'"} \
		~{"-samprate " + sample_ampl_rate} \
		~{"-lamprate " + locus_ampl_rate} \
		--ref_fasta ~{ref_genome_base} \
		~{"--amplicon_fasta " + ref_base} \
		--ampseq_export_format "~{ampseq_export_format}" \
		--poly_formula 'null' \
		--cigar_paths 'null' \
		--metadata metadata.csv \
		--join_by "Sample_id" \
		--Latitude "~{metadata_latitude_name}" \
		--Longitude "~{metadata_longitude_name}" \
		~{"--Variable1 \'" + metadata_variable1_name + "\'"} \
		~{"--Variable2 \'" + select_first([metadata_variable2_name, metadata_variable1_name]) + "\'"}
		
		echo 'Finished filtering ASVs!'

	>>>
	
	output {
		File markersTable_o = "references/markersTable.csv"
		File ampseq_object = "Results/~{out_prefix}.xlsx"
	}
	runtime {
		cpu: 1
		memory: "40 GiB"
		disks: "local-disk 10 HDD"
		bootDiskSizeGb: 10
		preemptible: 3 
		maxRetries: 1
		docker: 'jorgeamaya/ampseq:latest'
	}
}
