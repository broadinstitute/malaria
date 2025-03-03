version 1.0

task pca {
	input {
		Array[File] ampseq_excelfiles

		Float sample_ampl_rate = 0.2
		Float locus_ampl_rate = 0.2

		File ref_gff
		File ref_fasta
		File selected_checkboxes
		File? preserve_samples

		String join_by = "null"
		Boolean na_var_rm = true
		Boolean na_hap_rm = true

		Float ibd_thres = 0.9
		Boolean parallel = true
		Int ibd_ncol = 4
		String metadata_variable1_name = 'Country'
		String metadata_variable2_name = 'Run'
		String metadata_latitude_name = 'Latitude'
		String metadata_longitude_name = 'Longitude'

		Int nchunks = 100
		String hap_color_palette = "random"
		String poly_quantile = "null"
		String poly_formula = "null"
	}
	
	command <<<
		#FUTURE DEVELOPMENT
		# 1 REMOVE SAMPLE_METADATA DEPENDENCE. USER ALWAYS GET THAT FILE WRONG
		# 2 UPDATE THE SHINYAPP ACCORDINGLY
		# 3 DO NOT COPY THE FILES TO REFERENCE. USE NATIVE FILES.

		export TMPDIR=tmp
		set -euxo pipefail
		mkdir Reference
		mkdir Results
		mkdir Ampseq_Data 
		cp ~{ref_gff} Reference/.
		cp ~{ref_fasta} Reference/.
		cp ~{selected_checkboxes} Reference/.
		gsutil cp ~{sep=" " ampseq_excelfiles} Ampseq_Data/.
		
		#[TODO: Ask regarding placement of this statement]
		echo -e "\ngene_names_drug_resistance__\nPfDHFR\nPfMDR1\nPfDHPS\nPfKelch13C580Y\nPF3D7_1447900\ngene_ids_drug_resistance__\nPF3D7_0417200\nPF3D7_0523000\nPF3D7_0810800\nPF3D7_1343700\nPF3D7_1447900\ngene_names_diversity__\nCSP\nAMA1\nSERA2\nTRAP\ngene_ids_diversity__\nPF3D7_0304600\nPF3D7_1133400\nPF3D7_0207900\nPF3D7_1335900\n" >> ~{selected_checkboxes}

		Rscript /Code/MHap_Tertiary_Analysis_pipeline.R -fd /Code -ampseqe "Ampseq_Data" \
		-o "MHap_Profile" \
		-samprate ~{sample_ampl_rate} \
		-lamprate ~{locus_ampl_rate} \
		-gff ~{ref_gff} \
		-fasta ~{ref_fasta} \
		-join_by ~{join_by} \
		-preserve_samples ~{if defined(preserve_samples) then preserve_samples else "null"} \
		-Var1 ~{metadata_variable1_name} \
		-Var2 ~{metadata_variable2_name} \
		-Longitude ~{metadata_longitude_name} \
		-Latitude ~{metadata_latitude_name} \
		-na_var_rm ~{na_var_rm} \
		-na_hap_rm ~{na_hap_rm} \
		-ibd ~{ibd_thres} \
		-parallel ~{parallel} \
		-ibd_ncol ~{ibd_ncol} \
		-nchunks ~{nchunks} \
		-selected_checkboxes ~{selected_checkboxes} \
		-hap_color_palette "~{hap_color_palette}" \
		-poly_quantile ~{poly_quantile} \
		-pairwise_relatedness_table 'null' \
		-poly_formula "~{poly_formula}"	
	>>>

	output {
		File? ibs_connectivity_report = "Results/MHap_Profile_IBS_Connectivity_Report.html"
		File? ibd_connectivity_report = "Results/MHap_Profile_IBD_Connectivity_Report.html"
		#File? ibd_transmission_report = "Results/MHap_Profile_IBD_Transmission_Report.html"
		#File? coi_report = "Results/MHap_Profile_COI_Report.html"
	}

	runtime { 
		cpu: 16 
		memory: "40 GiB" 
		disks: "local-disk 10 HDD" 
		bootDiskSizeGb: 10 
		preemptible: 3 
		maxRetries: 1
		docker: 'jorgeamaya/pca:v_0_0_1'
	} 
}
