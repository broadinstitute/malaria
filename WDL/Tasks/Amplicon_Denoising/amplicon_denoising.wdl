version 1.0

task amplicon_denoising {
    input {
        Array[String] sample_ids
        Array[File] fastq1s
        Array[File] fastq2s
        Array[File] adaptor_rem1s
        Array[File] adaptor_rem2s
        Array[File] primer_rem1s
        Array[File] primer_rem2s
        File forward_primers_file
        File reverse_primers_file
        File reference_amplicons
        File? reference_amplicons_2
        Array[String] run_id

        #Variables for DADA2
        String Class = "parasite" # Name of the column in the metadata file for sample class
        String maxEE = "5,5" # Maximum expected error rate
        String trimRight = "0,0" # Number of bases to trim from the 3' end
        Int minLen = 30 # Minimum length of reads to retain
        String truncQ = "5,5" # Quality threshold for truncating reads
        String matchIDs = "0" # Whether to match IDs on fastqs
        Int max_consist = 10 # Maximum number of mismatches in overlap region
        Float omegaA = 0.000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001 # Alpha parameter for consensus quality score
        Int justConcatenate = 0 # Whether to just concatenate reads without merging
        Int maxMismatch = 0 # Maximum number of mismatches allowed during merging
        String saveRdata = "" # Whether to save the intermediate R data files

        #Variables for post-processing
        File? path_to_snv
        String strain = "3D7"
        String strain2 = "DD2"

        #Variables for ASV to CIGAR table conversion
        String min_reads = "0"
        String min_samples = "0"
        String max_snv_dist = "-1"
        String max_indel_dist = "-1"
        Boolean include_failed = false
        Boolean exclude_bimeras = false
        String polyN = "5"
    }

    command <<<
        # FUTURE DEVELOPMENT: 
        # 1 SOME FUNCTIONALITIES INSIDE RUNDADA2.R AND ASV_TO_CIGARTABLE.PY ARE NOT USED IN PIPELINE VERSION. CHECK IF THEY ARE NEEDED IN THE FUTURE.
        export TMPDIR=tmp
        set -euxo pipefail

        sample_ids_string=$(IFS=" "; echo "~{sep=' ' sample_ids}")
        fastq1s_string=$(IFS=" "; echo "~{sep=' ' fastq1s}")
        fastq2s_string=$(IFS=" "; echo "~{sep=' ' fastq2s}")
        adaptor_rem1s_string=$(IFS=" "; echo "~{sep=' ' adaptor_rem1s}")
        adaptor_rem2s_string=$(IFS=" "; echo "~{sep=' ' adaptor_rem2s}")
        primer_rem1s_string=$(IFS=" "; echo "~{sep=' ' primer_rem1s}")
        primer_rem2s_string=$(IFS=" "; echo "~{sep=' ' primer_rem2s}")

        mkdir -p Results

        ################################################################### 
        # Run DADA2 denoising                                             #
        ################################################################### 
        echo "Running DADA2 denoising..."

        Rscript /Code/runDADA2.R \
                    -b "${sample_ids_string}" \
                    -f1 "${fastq1s_string}" \
                    -f2 "${fastq2s_string}" \
                    -p1 "${primer_rem1s_string}" \
                    -p2 "${primer_rem2s_string}" \
                    -a1 "${adaptor_rem1s_string}" \
                    -a2 "${adaptor_rem2s_string}" \
                    -d "Results" \
                    -o "seqtab.tsv" \
                    -c "~{Class}" \
                    -ee "~{maxEE}" \
                    -tR "~{trimRight}" \
                    -mL "~{minLen}" \
                    -tQ "~{truncQ}" \
                    -id "~{matchIDs}" \
                    -mC "~{max_consist}" \
                    -wA "~{omegaA}" \
                    -jC "~{justConcatenate}" \
                    -mM "~{maxMismatch}" \
                    -s "~{saveRdata}" \
                    --bimera 

        echo "DADA2 denoising complete."

        ###################################################################
        # Run post-processing                                             #
        ###################################################################
        echo "Running post-processing..."

        Rscript /Code/postProc_dada2.R \
                    -s "Results/seqtab.tsv" \
                    -b "Results/ASVBimeras.txt" \
                    $(if [[ -f "~{path_to_snv}" ]]; then echo "-snv ~{path_to_snv}"; else echo "-snv No_File"; fi) \
                    --indel_filter "0.01" \
                    -o "Results" \
                    --fasta \
                    $(if [[ -f "~{reference_amplicons}" ]]; then echo "--reference_amplicons ~{reference_amplicons} --strain ~{strain}"; else echo "-no_ref"; fi) \
                    $(if [[ -f "~{reference_amplicons_2}" ]]; then echo "--reference_amplicons_2 ~{reference_amplicons_2} --strain2 ~{strain2}"; fi)

        echo "Post-processing complete."

        ###################################################################
        # Convert ASV to CIGAR tables                                     #
        ###################################################################
        echo "Converting ASV to CIGAR tables..."

        mkdir Results/Alignments
        python /Code/ASV_to_CIGARtable.py \
                    --reference_amplicons "~{reference_amplicons}" \
                    --seqtab "Results/seqtab.tsv" \
                    --ASVSeqs "Results/ASVSeqs.fasta" \
                    --ASVTable "Results/ASVTable.txt" \
                    --min_reads ~{min_reads} \
                    --min_samples ~{min_samples} \
                    --max_snv_dist ~{max_snv_dist} \
                    --max_indel_dist ~{max_indel_dist} \
                    --polyN ~{polyN} \
                    $(if ~{include_failed}; then echo "--include_failed"; fi) \
                    $(if ~{exclude_bimeras}; then echo "--exclude_bimeras"; fi) 

        echo "Conversion ASV to CIGAR tables complete."

    >>>

    output {
        #File seqtab_o = "Results/seqtab.tsv"
        File ASVBimeras_o = "Results/ASVBimeras.txt"
        File ASVTable_o = "Results/ASVTable.txt"
        File ASVSeqs_o = "Results/ASVSeqs.fasta"
        File CIGARVariants_Bfilter_o = "Results/CIGARVariants_Bfilter.out.tsv"
        File ASV_to_CIGAR_o = "Results/ASV_to_CIGAR.out.txt"
        File ZeroReadsSampleList_o = "Results/ZeroReadsSampleList.txt"

    }

    runtime {
        cpu: 4
        memory: "15 GiB"
        disks: "local-disk 10 HDD"
        bootDiskSizeGb: 10
        preemptible: 3
        maxRetries: 0
        docker: 'jorgeamaya/ampseq:latest'
    }
}