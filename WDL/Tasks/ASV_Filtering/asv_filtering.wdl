version 1.0

task asv_filtering {
	input {
		File markersTable
		File reference_genome
		File? panel_bedfile
		File? reference_amplicons
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
		File ASVSeqs
		File ASV_to_CIGAR
		File ZeroReadsSampleList
		File ReadAttrition

		# MHap ASV filtering thresholds 
		String sample_id_pat = '.'
		Int? min_abd = 10                                                                        # Minimum read depth to accept an ASV
		Float? min_ratio = 0.1                                                                   # Minimum ratio between major and minor alleles to accept a minor allele
		String? off_target_formula = "dVSITES_ij>=0.2"											 # Formula used to identify and remove off-target PCR products
		String? flanking_INDEL_formula = "flanking_INDEL==TRUE"						 # Formula to mask INDELs in the 3' and 5' ends of the ASVs
		Int? homopolymer_length = 5															     # Minimum number of single nucleotide tandem repeats to define homopolymeric regions in an ASV
		String? SNV_in_homopolymer_formula = "SNV_in_homopolymer==TRUE"              # Formula to mask SNVs in the homopolymeric regions of ASVs
		String? INDEL_in_homopolymer_formula = "INDEL_in_homopolymer==TRUE"          # Formula to mask INDELs in the homopolymeric regions of ASVs   
		String? bimera_formula = "bimera==TRUE\&h_ij>=0.66"                                       # Formula to detect bimeras
		String? PCR_errors_formula = "h_ij>=0.66\&h_ijminor>=0.66\&p_ij>=0.05"                   # Formula to detect other kinds of PCR errors
		

		Float? sample_ampl_rate = 0.3                                                            # Minimum proportion of amplified loci by a sample required to keep a sample
		Float? locus_ampl_rate = 0.3                                                             # Minimum proportion of amplified samples at a locus required to keep a locus
	}

	File ref_for_markers = select_first([panel_bedfile, reference_amplicons])
	###########################################
	# MHap - Define appropriate directories
	String wd = "Results/"
	String fd = "/Code/"
	String rd = "references/"
	String ref_genome_base = basename(reference_genome)
	String ref_base = if defined(reference_amplicons) then basename(select_first([reference_amplicons])) else ""
	# The directories below are subdirectories of "wd" (from MHap specs)
	String cigar_variants_dir = "cigar_variants/"
	String asv_table_dir = "asv_tables/"
	String asv2cigar_dir = "asv2cigar/"
	String asv_seq_dir = "asv_seq/"
	String zero_read_sample_list_dir = "zeroReadSampleList/"
        String read_attrition_dir = "read_attrition/"

	###########################################
	command <<<
		# FUTURE DEVELOPMENT
		# 1 REMOVE DEPENDENCE ON GSUTIL
		export TMPDIR=tmp
		set -euxo pipefail

		# Get amount of memory to use:
		mem_available=$(free -g | grep '^Mem' | awk '{print $2}')
		mem_max=$((mem_available-2))
		
		if [[ ${mem_max} -lt 2 ]] ; then
			mem_max=2
		fi

		echo "Available memory: ${mem_max} GiB"

		# Create directories - required for MHap script
		mkdir -p Results/
		mkdir -p Results/~{cigar_variants_dir}
		mkdir -p Results/~{asv_table_dir}
		mkdir -p Results/~{asv2cigar_dir}
		mkdir -p Results/~{asv_seq_dir}
		mkdir -p Results/~{zero_read_sample_list_dir}
		mkdir -p Results/~{read_attrition_dir}
		mkdir -p references/

		# Copy correct files to proper directories - required for MHap script
		gsutil cp ~{reference_amplicons} references/
		gsutil cp ~{reference_genome} references/
		~{"gsutil cp " + panel_bedfile + " references/"}
		gsutil cp ~{sample_metadata} Results/metadata.csv

		out_prefix=$(date +"%Y_%m_%d")
		gsutil cp ~{CIGARVariants} Results/~{cigar_variants_dir}/${out_prefix}_CIGARVariants_Bfilter.out.tsv
		gsutil cp ~{ASVTable} Results/~{asv_table_dir}
		gsutil cp ~{ASV_to_CIGAR} Results/~{asv2cigar_dir}
		gsutil cp ~{ASVSeqs} Results/~{asv_seq_dir}
		gsutil cp ~{ZeroReadsSampleList} Results/~{zero_read_sample_list_dir}
		gsutil cp ~{markersTable} references/markersTable.csv
		gsutil cp ~{ReadAttrition} Results/~{read_attrition_dir}/ #READATTRITION LOG HERE -> NEEDS INTEGRATION WITH XLSX FILE

		# Call MHap_Analysis_pipeline from Amplicon_TerraPipeline
		find . -type f
		#find /cromwell_root -type f
		#ls -l1 /cromwell_root
		root_dir=$(pwd)
		echo "Applying filters to ASVs..."
		Rscript /Code/MHap_Analysis_pipeline.R \
		-xmx "${mem_max}" \
		-fd "~{fd}" \
		-wd "$root_dir/~{wd}" \
		-rd "$root_dir/~{rd}" \
		-cigar_files '~{cigar_variants_dir}' \
		-asv_table_files '~{asv_table_dir}' \
		-asv2cigar_files '~{asv2cigar_dir}' \
		-asv_seq_files '~{asv_seq_dir}' \
		-zero_read_sample_list '~{zero_read_sample_list_dir}' \
		-o "${out_prefix}_ampseq_object_f" \
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
		File ampseq_object_o = glob("Results/*_ampseq_object_f.xlsx")[0]
	}
	runtime {
		cpu: 1
		memory: "100 GiB"
		disks: "local-disk 10 HDD"
		bootDiskSizeGb: 10
		preemptible: 3 
		maxRetries: 1
		docker: 'jorgeamaya/asvfilters:v_0_0_5'
	}
}
