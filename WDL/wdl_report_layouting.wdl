version 1.0

workflow report_layouting {
	input {
		Array[File] cigar_files
		Array[File] metadata_files
		
		Int nTasks = 50

		String cigar_paths = "NaN"
		String cigar_dir = "cigar_dir"
		String ampseq_jsonfile = "NaN"
		String ampseq_excelfile = "NaN"

		String sample_id_pattern = "^ID"
		File markers 
		Int min_abd = 10
		Float min_ratio = 0.1

		Boolean PerformanceReport = false

		Float sample_ampl_rate = 0.75
		Float locus_ampl_rate = 0.75

		Boolean Drug_Surveillance_Report = true
		Boolean Variants_of_Interest_Report = false

		File ref_gff
		File ref_fasta
		File reference_alleles
		File selected_checkboxes
	
		String join_by = "Sample_id"
		Boolean na_var_rm = true
		Boolean na_hap_rm = true
		String drugs = "Artemisinin,Chloroquine,Pyrimethamine,Sulfadoxine,Lumefantrine,Mefloquine"
		Boolean include_all_drug_markers = true

		String ibd_thres = "NaN"
		Boolean parallel = true
		Int ibd_ncol = 4
		String? pop_levels = "null"

		Int nchunks = 500
	}
	
	call report_layouting_process {
		input:
			cigar_files = cigar_files,
			metadata_files = metadata_files,
			nTasks = nTasks,
			cigar_paths = cigar_paths,
			cigar_dir = cigar_dir,
			ampseq_jsonfile = ampseq_jsonfile,
			ampseq_excelfile = ampseq_excelfile,
			sample_id_pattern = sample_id_pattern,
			markers = markers,
			min_abd = min_abd,
			min_ratio = min_ratio,
			PerformanceReport = PerformanceReport,
			sample_ampl_rate = sample_ampl_rate,
			locus_ampl_rate = locus_ampl_rate,
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
			nchunks = nchunks
	}

	output {
		File? html_report_f = report_layouting_process.html_report
	}
}

task report_layouting_process {
	input {
		Array[File] cigar_files
		Array[File] metadata_files

		Int nTasks

		String cigar_paths
		String cigar_dir
		String ampseq_jsonfile
		String ampseq_excelfile

		String sample_id_pattern
		File markers 
		Int min_abd
		Float min_ratio

		Boolean PerformanceReport

		Float sample_ampl_rate
		Float locus_ampl_rate

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

		Int nchunks
	}
	
	command <<<
		set -euxo pipefail
		mkdir cigar_dir
		mkdir Reference
		mkdir Results
		cp ~{sep = ' ' cigar_files} cigar_dir
		#gsutil -m cp -r ~{sep = ' ' cigar_files} cigar_dir/
		cp ~{ref_gff} Reference/.
		cp ~{ref_fasta} Reference/.
		cp ~{reference_alleles} Reference/.
		cp ~{markers} Reference/.
		echo "CIGAR TABLES"
		ls cigar_dir
		
		echo "Sample_id,Geo_Level,Temp_Level,Longitude,Latitude" > metadata.csv
		cat ~{sep = ' ' metadata_files} >> metadata.csv

		Rscript /Code/MHap_Analysis_pipeline.R -fd /Code -cigar_paths ~{cigar_paths} \
		-cigar_files "cigar_dir" \
		-ampseqj ~{ampseq_jsonfile} \
		-ampseqe ~{ampseq_excelfile} \
		-o "MHap_Profile" \
		-sample_id_pattern ~{sample_id_pattern} \
		-markers ~{markers} \
		-min_abd ~{min_abd} \
		-min_ratio ~{min_ratio} \
		-samprate ~{sample_ampl_rate} \
		-lamprate ~{locus_ampl_rate} \
		-PerformanceReport ~{PerformanceReport} \
		-Drug_Surveillance_Report ~{Drug_Surveillance_Report} \
		-Variants_of_Interest_Report ~{Variants_of_Interest_Report} \
		-gff ~{ref_gff} \
		-fasta ~{ref_fasta} \
		-reference_alleles ~{reference_alleles} \
		-metadata metadata.csv \
		-join_by ~{join_by} \
		-Var1 Geo_Level \
		-Var2 Temp_Level \
		-Longitude Longitude \
		-Latitude Latitude \
		-na_var_rm ~{na_var_rm} \
		-na_hap_rm ~{na_hap_rm} \
		-drugs ~{drugs} \
		-include_all_drug_markers ~{include_all_drug_markers} \
		-t ~{nTasks} \
		-tid 1 \
		-ibd ~{ibd_thres} \
		-parallel ~{parallel} \
		-ibd_ncol ~{ibd_ncol} \
		-pop_levels ~{pop_levels} \
		-nchunks ~{nchunks} \
		-selected_checkboxes ~{selected_checkboxes}
	
		ls 
		ls Results/
	>>>

	output {
		File? html_report = "Results/MHap_Profile_DRS_Report.html"
	}

	runtime { 
		cpu: 1 
		memory: "15 GiB" 
		disks: "local-disk 10 HDD" 
		bootDiskSizeGb: 10 
		preemptible: 3 
		maxRetries: 1 
		docker: 'jorgeamaya/mhap' 
	} 
}
