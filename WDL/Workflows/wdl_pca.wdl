version 1.0

import "../Tasks/PCA/pca_process.wdl" as pca_process_t

workflow plot_pca {
	input {
		Array[File] ampseq_excelfiles

		File ref_gff
		File ref_fasta
		File selected_checkboxes
		File? preserve_samples
		File? reference_samples
	}
	
	call pca_process_t.pca as t_001_pca_process {
		input:
			ampseq_excelfiles = ampseq_excelfiles,
			ref_gff = ref_gff,
			ref_fasta = ref_fasta,
			selected_checkboxes = selected_checkboxes,
			preserve_samples = preserve_samples,
			reference_samples = reference_samples

	}

	output {
		File? ibs_connectivity_report_f = t_001_pca_process.ibs_connectivity_report
		File? ibd_connectivity_report_f = t_001_pca_process.ibd_connectivity_report
		#File? ibd_transmission_report_f = t_001_pca_process.ibd_transmission_report
		#File? coi_report_f = t_001_report_layouting_process.coi_report
	}
}

