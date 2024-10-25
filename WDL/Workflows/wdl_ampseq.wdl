version 1.0

import "../Tasks/Prepare_Reference_Files/prepare_reference_files.wdl" as prepare_reference_files_t
import "../Tasks/Cutadapters/cutadapters.wdl" as cutadapters_t
import "../Tasks/Trimprimers/trimprimers.wdl" as trimprimers_t
import "../Tasks/Amplicon_Denoising/amplicon_denoising.wdl" as amplicon_denoising_t
import "../Tasks/ASV_Filtering/asv_filtering.wdl" as asv_filtering_t

#import "../Tasks/amplicon_demultiplexing.wdl" as amplicon_demultiplexing_t
#import "../Tasks/contamination_detection.wdl" as contamination_detection_t

workflow ampseq {
    input {    
        # Required input files
        Array[File] fastq1s
        Array[File] fastq2s
        Array[String] sample_ids
        File forward_primers_file
        File reverse_primers_file
        File reference_genome
        Array[String] run_id

        # Optional reference files
        File? reference_amplicons
        File? reference_amplicons_2 
        File? panel_bedfile
        File? markersTable
        File? path_to_snv

        String pattern_fw = "_L001_R1_001.fastq.gz"
        String pattern_rv = "_L001_R2_001.fastq.gz"
        
        # Optional file for decontamination pipeline
        File? barcodes_index

        # Parameters for adaptor and primer removal
        Int trim_galore_qvalue = 5
        Int trim_galore_length = 20

        # Commands for the contamination detection pipeline
        Int minreads_threshold = 1000
        Float contamination_threshold = 0.5
        String verbose = "False"
        Boolean run_demultiplexing = false

        #TMP FOR TESTING OF CONTAMINATION DETECTION PIPELINE
        #File? path_to_flist
        #File? config_json
    }

    call prepare_reference_files_t.prepare_reference_files as t_001_prepare_reference_files {
        input:
            sample_ids = sample_ids,
            forward_primers_file = forward_primers_file,
            reverse_primers_file = reverse_primers_file,
            reference_amplicons = reference_amplicons,
            reference_amplicons_2 = reference_amplicons_2,
            path_to_snv = path_to_snv,
            reference_genome = reference_genome,
            panel_bedfile = panel_bedfile,
            pattern_fw = pattern_fw,
            pattern_rv = pattern_rv,
            minreads_threshold = minreads_threshold,
            contamination_threshold = contamination_threshold,
            verbose = verbose,
            barcodes_index = barcodes_index
    }

#    if(defined(barcodes_index)) {
#        call contamination_detection_t.contamination_detection {
#            input: 
#                #TMP FOR TESTING OF CONTAMINATION DETECTION PIPELINE
#                config_json = t_001_prepare_reference_files.config_json_out,
#                path_to_flist = t_001_prepare_reference_files.path_to_flist_o,
#                fastq1s = fastq1s,
#                fastq2s = fastq2s,
#                forward_primers_file = forward_primers_file,
#                reverse_primers_file = reverse_primers_file,
#                barcodes_index = barcodes_index
#        }
#    }
#
    #if(!run_demultiplexing) {

    scatter (indx1 in range(length(fastq1s))) {
        call cutadapters_t.cutadapters as t_002_cutadapters {
            input:
                fastq1 = fastq1s[indx1],
                fastq2 = fastq2s[indx1],
                trim_galore_qvalue = trim_galore_qvalue,
                trim_galore_length = trim_galore_length,
        }

        call trimprimers_t.trimprimers as t_003_trimprimers {
            input:
                trimmed_fastq1 = t_002_cutadapters.fastq1_noadapters,
                trimmed_fastq2 = t_002_cutadapters.fastq2_noadapters,
                forward_primers_file = forward_primers_file,
                reverse_primers_file = reverse_primers_file,
        }

    }


#    if(run_demultiplexing) {
#        call ampseq_pipeline_demult {
#            input: 
#                config_json = t_001_prepare_reference_files.config_json_out,
#                path_to_flist = path_to_flist,
#                fastq1s = fatsq1s,
#                fastq2s = fatsq2s,
#                forward_primers_file = forward_primers_file,
#                reverse_primers_file = reverse_primers_file,
#                reference = t_001_prepare_reference_files.reference_out,
#                reference2 = reference_amplicons_2,
#                path_to_snv = path_to_snv
#        }
#    }

    call amplicon_denoising_t.amplicon_denoising as t_004_amplicon_denoising {
        input:
            fastq1s = fastq1s,
            fastq2s = fastq2s,
            forward_primers_file = forward_primers_file,
            reverse_primers_file = reverse_primers_file,
            reference = t_001_prepare_reference_files.reference_out,
            reference2 = reference_amplicons_2,
            run_id = run_id,
            path_to_snv = path_to_snv,
            primer_rem = flatten([t_003_trimprimers.fastq1_no_primers, t_003_trimprimers.fastq2_no_primers]),
            adaptor_rem = flatten([t_002_cutadapters.fastq1_noadapters, t_002_cutadapters.fastq2_noadapters])
    }
#
#    call asv_filtering_t.asv_filtering {
#        input: 
#            reference = reference_amplicons,
#            reference_genome = reference_genome,
#            panel_bedfile = panel_bedfile,
#            markersTable = markersTable,        
#            CIGARVariants = t_004_amplicon_denoising.CIGARVariants_Bfilter,
#            ASVTable = t_004_amplicon_denoising.ASVTable,
#            ASVSeqs = t_004_amplicon_denoising.ASVSeqs,
#            ASV_to_CIGAR = t_004_amplicon_denoising.ASV_to_CIGAR,
#            ZeroReadsSampleList = t_004_amplicon_denoising.ZeroReadsSampleList
#    }

    output {
        # ASV Filtering
        #File ampseq_object_f = asv_filtering.ampseq_object

        # Keep this variables for testing purposes
        #File panel_reference_fasta_f = t_001_prepare_reference_files.reference_out
        
        # DADA2
        # File ASVBimeras_f = ampseq_pipeline_denoise.ASVBimeras
        # File seqtab_f = ampseq_pipeline_denoise.seqtab

        # PostProc_DADA2
        # File ASVTable_f = ampseq_pipeline_denoise.ASVTable
        # File ASVSeqs_f = ampseq_pipeline_denoise.ASVSeqs

        # ASV_to_CIGAR
        # File CIGARVariants_Bfilter_f = ampseq_pipeline_denoise.CIGARVariants_Bfilter
        # File ASV_to_CIGAR_f = ampseq_pipeline_denoise.ASV_to_CIGAR
        # File ZeroReadsSampleList_f = ampseq_pipeline_denoise.ZeroReadsSampleList

        # File markersTable_f = ampseq_pipeline_asv_filtering.markersTable_o

        ###REMOVE THIS VARIABLES AFTER TESTING###
        #File path_to_flist_o_f = t_001_prepare_reference_files.path_to_flist_o
        #File? contamination_detection_missing_files_f = contamination_detection.missing_files
        #File? contamination_detection_sample_cards_f = contamination_detection.contamination_detection_sample_cards
        # File? contamination_detection_report_f = contamination_detection.contamination_detection_report

        Array[File] cutadaptersout_f_f = t_002_cutadapters.fastq1_noadapters
        Array[File] cutadaptersout_r_f = t_002_cutadapters.fastq2_noadapters
        Array[File] trimprimersout_f_f = t_003_trimprimers.fastq1_no_primers
        Array[File] trimprimersout_r_f = t_003_trimprimers.fastq2_no_primers
    }
}

