version 1.0

task report_layouting_process {
	input {
		File ampseq_excelfile

		Boolean PerformanceReport

		Float sample_ampl_rate = 0.2
		Float locus_ampl_rate = 0.2

		Boolean Drug_Surveillance_Report
		Boolean Variants_of_Interest_Report

		File ref_gff
		File ref_fasta
		File reference_alleles
		File selected_checkboxes
	
		String join_by
		Boolean na_var_rm
		Boolean na_hap_rm
		String drugs
		Boolean include_all_drug_markers

		String ibd_thres
		Boolean parallel
		Int ibd_ncol
		String? pop_levels
		String metadata_variable1_name = 'Country'
		String metadata_variable2_name = 'Run'
		String metadata_latitude_name = 'Latitude'
		String metadata_longitude_name = 'Longitude'

		Int nchunks
		String hap_color_palette
		String poly_quantile
		String poly_formula
	}
	
	command <<<
		#FUTURE DEVELOPMENT
		# 1 REMOVE SAMPLE_METADATA DEPENDENCE. USER ALWAYS GET THAT FILE WRONG
		# 2 UPDATE THE SHINYAPP ACCORDINGLY
		# 3 DO NOT COPY THE FILES TO REFERENCE. USE NATIVE FILES.

		set -euxo pipefail
		mkdir Reference
		mkdir Results
		cp ~{ref_gff} Reference/.
		cp ~{ref_fasta} Reference/.
		cp ~{reference_alleles} Reference/.
		cp ~{selected_checkboxes} Reference/.
		echo "CIGAR TABLES"
		
		#[TODO: Ask regarding placement of this statement]
		echo -e "\ngene_names_drug_resistance__\nPfDHFR\nPfMDR1\nPfDHPS\nPfKelch13C580Y\nPF3D7_1447900\ngene_ids_drug_resistance__\nPF3D7_0417200\nPF3D7_0523000\nPF3D7_0810800\nPF3D7_1343700\nPF3D7_1447900\ngene_names_diversity__\nCSP\nAMA1\nSERA2\nTRAP\ngene_ids_diversity__\nPF3D7_0304600\nPF3D7_1133400\nPF3D7_0207900\nPF3D7_1335900\n" >> ~{selected_checkboxes}

		Rscript /Code/MHap_Tertiary_Analysis_pipeline.R -fd /Code -ampseqe ~{ampseq_excelfile} \
		-o "MHap_Profile" \
		-samprate ~{sample_ampl_rate} \
		-lamprate ~{locus_ampl_rate} \
		-PerformanceReport ~{PerformanceReport} \
		-Drug_Surveillance_Report ~{Drug_Surveillance_Report} \
		-Variants_of_Interest_Report ~{Variants_of_Interest_Report} \
		-gff ~{ref_gff} \
		-fasta ~{ref_fasta} \
		-reference_alleles ~{reference_alleles} \
		-join_by ~{join_by} \
		-Var1 ~{metadata_variable1_name} \
		-Var2 ~{metadata_variable2_name} \
		-Longitude ~{metadata_longitude_name} \
		-Latitude ~{metadata_latitude_name} \
		-na_var_rm ~{na_var_rm} \
		-na_hap_rm ~{na_hap_rm} \
		-drugs ~{drugs} \
		-include_all_drug_markers ~{include_all_drug_markers} \
		-ibd ~{ibd_thres} \
		-parallel ~{parallel} \
		-ibd_ncol ~{ibd_ncol} \
		-pop_levels ~{pop_levels} \
		-nchunks ~{nchunks} \
		-selected_checkboxes ~{selected_checkboxes} \
		-hap_color_palette "~{hap_color_palette}" \
		-poly_quantile ~{poly_quantile} \
		-pairwise_relatedness_table 'null' \
		-poly_formula "~{poly_formula}"
	
	>>>

	output {
		#File? intermediate_report = "Results/metadata_intermediate.csv"
		#File? drs_report = "Results/MHap_Profile_DRS_Report.html"
		File? drs_minimal_report = "Results/MHap_Profile_DRS_Minimal_Report.html"
		#File? coi_report = "Results/MHap_Profile_COI_Report.html"
		#File? ibd_connectivity_report = "Results/MHap_Profile_IBD_Connectivity_Report.html"
		#File? ibd_transmssion_report = "Results/MHap_Profile_IBD_Transmission_Report.html"
		#File? performance_report = "Results/MHap_Profile_Performance_Report.html"
	}

	runtime { 
		cpu: 16 
		memory: "40 GiB" 
		disks: "local-disk 10 HDD" 
		bootDiskSizeGb: 10 
		preemptible: 3 
		maxRetries: 1
		docker: 'jorgeamaya/mhap:v_0_0_2'
	} 
}
