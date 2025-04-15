version 1.0

import "../Tasks/Prepare_Reference_Files/prepare_reference_files.wdl" as prepare_reference_files_t
import "../Tasks/File_Validation/validate_inputs.wdl" as validate_inputs_t
import "../Tasks/Contamination_Detection/contamination_detection.wdl" as contamination_detection_t
import "../Tasks/Cutadapters/cutadapters.wdl" as cutadapters_t
import "../Tasks/Trimprimers/trimprimers.wdl" as trimprimers_t
import "../Tasks/Amplicon_Denoising/amplicon_denoising.wdl" as amplicon_denoising_t
import "../Tasks/ASV_Filtering/asv_filtering.wdl" as asv_filtering_t
#import "../Tasks/amplicon_demultiplexing.wdl" as amplicon_demultiplexing_t

workflow ampseq {
    input {    
        # Required input files
        Array[File] fastq1s
        Array[File] fastq2s
        Array[String]? sample_ids
        File reference_genome
        File panel_info

        # Optional reference files
        File? forward_primers_file
        File? reverse_primers_file
        File? reference_amplicons
        File? reference_amplicons_2 
        File? markers_table
        File? path_to_snv

        # Optional file for contamination detection pipeline
        File? barcodes_index

        # Parameters for Denoising and ASV filtering
        File sample_metadata
    }
    #FUTURE DEVELOPMENT
    # 1 CLEAN DOCKERFILES
    # 2 REMOVE File sample_metadata. USERS WILL ALWAYS GET THIS WRONG.
    # 3 ADD FEATURE TO DOWN SAMPLE FILES
    # 4 DISCONTINUE THE MARKERS TABLE?

    call prepare_reference_files_t.prepare_reference_files as t_000_prepare_reference_files {
        input:
            panel_info = panel_info,
            reference_genome = reference_genome,
            reference_amplicons = reference_amplicons,
            reference_amplicons_2 = reference_amplicons_2,
            forward_primers_file = forward_primers_file,
            reverse_primers_file = reverse_primers_file,
            path_to_snv = path_to_snv
    }

    scatter (indx1 in range(length(fastq1s))) {
        call validate_inputs_t.validate_fastqs as t_001_validate_fastqs {
            input:
                fastq1 = fastq1s[indx1],
                fastq2 = fastq2s[indx1]
        }

        call cutadapters_t.cutadapters as t_002_cutadapters {
            input:
                fastq1 = t_001_validate_fastqs.fastq1_o,
                fastq2 = t_001_validate_fastqs.fastq2_o
        }

        call trimprimers_t.trimprimers as t_003_trimprimers {
            input:
                trimmed_fastq1 = t_002_cutadapters.fastq1_noadapters_o,
                trimmed_fastq2 = t_002_cutadapters.fastq2_noadapters_o,
                forward_primers_file = select_first([forward_primers_file, t_000_prepare_reference_files.forward_primers_o]),
                reverse_primers_file = select_first([reverse_primers_file, t_000_prepare_reference_files.reverse_primers_o])
        }
    }

    if(defined(barcodes_index)) {
        call contamination_detection_t.contamination_detection as t_0001_contamination_detection {
            input: 
                adaptor_rem1s = t_002_cutadapters.fastq1_noadapters_o,
                adaptor_rem2s = t_002_cutadapters.fastq2_noadapters_o,
                forward_primers_file = t_000_prepare_reference_files.forward_primers_CI_o,
                reverse_primers_file = t_000_prepare_reference_files.reverse_primers_CI_o,
                barcodes_index = barcodes_index
        }
    }

#    if(run_demultiplexing) {
#        call ampseq_pipeline_demult {
#            input: 
#                config_json = t_000_prepare_reference_files.config_json_out,
#                path_to_flist = path_to_flist,
#                fastq1s = fatsq1s,
#                fastq2s = fatsq2s,
#                forward_primers_file = forward_primers_file,
#                reverse_primers_file = reverse_primers_file,
#                reference = t_000_prepare_reference_files.reference_out,
#                reference2 = reference_amplicons_2,
#                path_to_snv = path_to_snv
#        }
#    }

#    call amplicon_denoising_t.amplicon_denoising as t_004_amplicon_denoising {
#        input:
#            sample_ids = sample_ids,
#            fastq1s = fastq1s,
#            fastq2s = fastq2s,
#            adaptor_rem1s = t_002_cutadapters.fastq1_noadapters_o,
#            adaptor_rem2s = t_002_cutadapters.fastq2_noadapters_o,
#            primer_rem1s = t_003_trimprimers.fastq1_no_primers_o,
#            primer_rem2s = t_003_trimprimers.fastq2_no_primers_o,
#            forward_primers_file = select_first([forward_primers_file, t_000_prepare_reference_files.forward_primers_o]),
#            reverse_primers_file = select_first([reverse_primers_file, t_000_prepare_reference_files.reverse_primers_o]),
#            reference_amplicons = t_000_prepare_reference_files.reference_o,
#            reference_amplicons_2 = reference_amplicons_2,
#            path_to_snv = path_to_snv
#    }
#
#    call asv_filtering_t.asv_filtering as t_005_asv_filtering {
#        input:
#            sample_metadata = sample_metadata,
#            panel_bedfile = t_000_prepare_reference_files.panel_bedfile_o,
#            reference_amplicons = t_000_prepare_reference_files.reference_o,
#            markersTable = t_000_prepare_reference_files.markers_table_o,      
#            reference_genome = reference_genome,
#            CIGARVariants = t_004_amplicon_denoising.CIGARVariants_Bfilter_o,
#            ASVTable = t_004_amplicon_denoising.ASVTable_o,
#            ASVSeqs = t_004_amplicon_denoising.ASVSeqs_o,
#            ASV_to_CIGAR = t_004_amplicon_denoising.ASV_to_CIGAR_o,
#            ZeroReadsSampleList = t_004_amplicon_denoising.ZeroReadsSampleList_o,
#            ReadAttrition = t_004_amplicon_denoising.ReadAttrition_o
#    }
#
    output {
        # PREPARE REFERENCE FILES
        #File reference_out = t_000_prepare_reference_files.reference_o
        #File panel_bedfile_out = t_000_prepare_reference_files.panel_bedfile_o
        #File forward_primers_out = t_000_prepare_reference_files.forward_primers_o
        #File reverse_primers_out = t_000_prepare_reference_files.reverse_primers_o
        #File markers_table_out = t_000_prepare_reference_files.markers_table_o

        # VALIDATE FASTQS
        #Array[File] fastq1_out = t_001_validate_fastqs.fastq1_o
        #Array[File] fastq2_out = t_001_validate_fastqs.fastq2_o
        
        # CUTADAPTERS AND TRIMPRIMERS
        #Array[File] cutadaptersout_f_out = t_002_cutadapters.fastq1_noadapters_o
        #Array[File] cutadaptersout_r_out = t_002_cutadapters.fastq2_noadapters_o
        #Array[File] trimprimersout_f_out = t_003_trimprimers.fastq1_no_primers_o
        #Array[File] trimprimersout_r_out = t_003_trimprimers.fastq2_no_primers_o

        # CONTAMINATION DETECTION
        File? contamination_detection_missing_files_f = t_0001_contamination_detection.missing_files
        File? contamination_detection_sample_cards_f = t_0001_contamination_detection.contamination_detection_sample_cards
        File? contamination_detection_report_f = t_0001_contamination_detection.contamination_detection_report

        # DADA2
        #File ASVBimeras_f = t_004_amplicon_denoising.ASVBimeras_o
        #File ReadAttrition_f = t_004_amplicon_denoising.ReadAttrition_o

        # PostProc_DADA2
        #File ASVTable_f = t_004_amplicon_denoising.ASVTable_o
        #File ASVSeqs_f = t_004_amplicon_denoising.ASVSeqs_o

        # ASV_to_CIGAR
        #File CIGARVariants_Bfilter_f = t_004_amplicon_denoising.CIGARVariants_Bfilter_o
        #File ASV_to_CIGAR_f = t_004_amplicon_denoising.ASV_to_CIGAR_o
        #File ZeroReadsSampleList_f = t_004_amplicon_denoising.ZeroReadsSampleList_o

        # ASV Filtering
#        File ampseq_object_f = t_005_asv_filtering.ampseq_object_o
    }
}

