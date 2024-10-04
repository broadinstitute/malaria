version 1.0

task amplicon_denoising {
	input {
		Array[File] fastq1s
		Array[File] fastq2s
		File path_to_flist
		File forward_primers_file
		File reverse_primers_file
		File reference
		File? reference2
		File? path_to_snv
		Array[File]? primer_rem
		Array[File]? adaptor_rem
		Array[String] run_id
		File config_json
	}

	Map[String, String] configs = read_json(config_json)
	String pattern_fw = configs["pattern_fw"]
	String pattern_rv = configs["pattern_rv"]

	command <<<
	set -euxo pipefail
	cat ~{config_json}

	###################################################################
	##Copy files to the working directory and run the AmpSeq pipeline##
	###################################################################
	mkdir -p Results/
	mkdir -p Results/PrimerRem
	mkdir -p Results/AdaptorRem
	mkdir -p Results/PostProc_DADA2
	mkdir -p Results/ASV_to_CIGAR
	mkdir -p Results/alignments

	mkdir fq_dir
	mkdir references
	gsutil -m cp -r ~{sep = ' ' fastq1s} fq_dir/

	# Match suggested pattern_fw baked into run_DADA2
	# [TODO: Fix reliance on specific suffixes]
	for file in fq_dir/*~{pattern_fw}; do
		mv -- "$file" "${file%~{pattern_fw}}_L001_R1_001.fastq.gz"
	done

	# Match suggested pattern_rv baked into run_DADA2
	# [TODO: Fix reliance on specific suffixes]
	gsutil -m cp -r ~{sep = ' ' fastq2s} fq_dir/
	for file in fq_dir/*~{pattern_rv}; do
		mv -- "$file" "${file%~{pattern_rv}}_L001_R2_001.fastq.gz"
	done
	
	gsutil cp ~{path_to_flist} references/samples.txt
	gsutil -m cp -r ~{forward_primers_file} references/
	gsutil -m cp -r ~{reverse_primers_file} references/
	gsutil -m cp -r ~{sep = ' ' primer_rem} Results/PrimerRem/
	gsutil -m cp -r ~{sep = ' ' adaptor_rem} Results/AdaptorRem/
	gsutil cp ~{reference} references/
	gsutil cp ~{path_to_flist} references/samples.txt

	~{"gsutil cp " + reference2 + " references/reference2.fasta"}
	~{"gsutil cp " + path_to_snv + " references/snv_filter.tsv"}

	if [ -e "Results/PrimerRem/mixed_nop_prim_meta.tsv" ]; then
		echo "Demultiplex performed in the data. Some read pairs assumed to be too short to overlap and be merged."
		python /Code/Amplicon_TerraPipeline.py --config ~{config_json} --terra --dada2 --demultiplexed
	else
		echo "Demultiplex not performed in the data. Read pairs assumed to be long enough to overlap and be merged."
		python /Code/Amplicon_TerraPipeline.py --config ~{config_json} --terra --dada2 
	fi

	# Run postproc_DADA2
	echo "Running post processing of DADA2 results..."
	python /Code/Amplicon_TerraPipeline.py --config ~{config_json} --terra --postproc_dada2

	# Run ASV_to_CIGAR
	echo "Converting ASV to CIGAR tables..."
	find . -type f
	python /Code/Amplicon_TerraPipeline.py --config ~{config_json} --terra --asv_to_cigar
	find . -type f
	>>>

	output {
		File ASVBimeras = "Results/ASVBimeras.txt"
		File CIGARVariants_Bfilter = "Results/CIGARVariants_Bfilter.out.tsv"
		File ASV_to_CIGAR = "Results/ASV_to_CIGAR/ASV_to_CIGAR.out.txt"
		File ZeroReadsSampleList = "Results/ASV_to_CIGAR/ZeroReadsSampleList.txt"

		File seqtab = "Results/seqtab.tsv"
		File ASVTable = "Results/PostProc_DADA2/ASVTable.txt"
		File ASVSeqs = "Results/PostProc_DADA2/ASVSeqs.fasta"
	}

	runtime {
		cpu: 1
		memory: "15 GiB"
		disks: "local-disk 10 HDD"
		bootDiskSizeGb: 10
		preemptible: 0
		maxRetries: 1
		docker: 'jorgeamaya/ampseq:v1.0'
	}
}
