version 1.0

import "../Tasks/Report_Layouting/report_layouting_process.wdl" as report_layouting_process_t

workflow report_layouting {
	input {
		File ampseq_excelfile

		Boolean PerformanceReport = true
		Boolean Drug_Surveillance_Report = true
		Boolean Variants_of_Interest_Report = false

		File ref_gff
		File ref_fasta
		File reference_alleles
		File selected_checkboxes
	
		String join_by = "null"
		Boolean na_var_rm = true
		Boolean na_hap_rm = true
		String drugs = "Artemisinin,Chloroquine,Pyrimethamine,Sulfadoxine,Lumefantrine,Mefloquine"
		Boolean include_all_drug_markers = false

		String ibd_thres = "null"
		Boolean parallel = true
		Int ibd_ncol = 4
		String? pop_levels = "null"

		String metadata_variable1_name = 'Country'
		String metadata_variable2_name = 'Run'
		String metadata_longitude_name = 'Longitude'
		String metadata_latitude_name = 'Latitude'

		Int nchunks = 100
		String hap_color_palette = "random"
		String poly_quantile = "null"
		String poly_formula = "null"

	}
	
	call report_layouting_process_t.report_layouting_process as t_001_report_layouting_process {
		input:
			ampseq_excelfile = ampseq_excelfile,
			PerformanceReport = PerformanceReport,
			Drug_Surveillance_Report = Drug_Surveillance_Report,
			Variants_of_Interest_Report = Variants_of_Interest_Report,
			ref_gff = ref_gff,
			ref_fasta = ref_fasta,
			reference_alleles = reference_alleles,
			selected_checkboxes = selected_checkboxes,
			join_by = join_by,
			na_var_rm = na_var_rm,
			na_hap_rm = na_hap_rm,
			drugs = drugs,
			include_all_drug_markers = include_all_drug_markers,
			ibd_thres = ibd_thres,
			parallel = parallel,
			ibd_ncol = ibd_ncol,
			pop_levels = pop_levels,
			nchunks = nchunks,
			hap_color_palette = hap_color_palette,
			poly_quantile = poly_quantile,
			poly_formula = poly_formula,
			metadata_longitude_name = metadata_longitude_name,
			metadata_latitude_name = metadata_latitude_name,
			metadata_variable1_name = metadata_variable1_name,
			metadata_variable2_name = metadata_variable2_name
	}

	output {
		#File? intermediate_report_f = t_001_report_layouting_process.intermediate_report
		#File? drs_report_f = t_001_report_layouting_process.drs_report
		File? drs_minimal_report_f = t_001_report_layouting_process.drs_minimal_report
		#File? coi_report_f = t_001_report_layouting_process.coi_report
		#File? ibd_connectivity_report_f = t_001_report_layouting_process.ibd_connectivity_report
		#File? ibd_transmssion_report_f = t_001_report_layouting_process.ibd_transmssion_report
		#File? performance_report_f = t_001_report_layouting_process.performance_report
	}
}

