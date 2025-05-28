version 1.0

import "../Tasks/Report_Layouting/report_layouting_process.wdl" as report_layouting_process_t

workflow report_layouting {
	input {
		File ampseq_excelfile
		File ref_gff
		File ref_fasta
		File reference_alleles
		File selected_checkboxes
	}
	
	call report_layouting_process_t.report_layouting_process as t_001_report_layouting_process {
		input:
			ampseq_excelfile = ampseq_excelfile,
			ref_gff = ref_gff,
			ref_fasta = ref_fasta,
			reference_alleles = reference_alleles,
			selected_checkboxes = selected_checkboxes
	}

	output {
		#File? intermediate_report_f = t_001_report_layouting_process.intermediate_report
		File? drs_report_f = t_001_report_layouting_process.drs_report
		File? drs_minimal_report_f = t_001_report_layouting_process.drs_minimal_report
		#File? coi_report_f = t_001_report_layouting_process.coi_report
		#File? ibd_connectivity_report_f = t_001_report_layouting_process.ibd_connectivity_report
		#File? ibd_transmssion_report_f = t_001_report_layouting_process.ibd_transmssion_report
		#File? performance_report_f = t_001_report_layouting_process.performance_report
	}
}

