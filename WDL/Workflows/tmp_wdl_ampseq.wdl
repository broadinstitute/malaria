version 1.0

#import "../Tasks/Prepare_Reference_Files/prepare_reference_files.wdl" as prepare_reference_files_t
#import "../Tasks/Cutadapters/cutadapters.wdl" as cutadapters_t
#import "../Tasks/Trimprimers/trimprimers.wdl" as trimprimers_t
#import "../Tasks/Amplicon_Denoising/amplicon_denoising.wdl" as amplicon_denoising_t
import "../Tasks/ASV_Filtering/asv_filtering.wdl" as asv_filtering_t

#import "../Tasks/amplicon_demultiplexing.wdl" as amplicon_demultiplexing_t
#import "../Tasks/contamination_detection.wdl" as contamination_detection_t

workflow ampseq {
    input {    
        # Required input files
        Array[File] fastq1s
        Array[File] fastq2s
        Array[String] sample_ids
        File reference_genome
        Array[String] run_id
        File panel_info

        # Optional reference files
        File? forward_primers_file
        File? reverse_primers_file
        File? reference_amplicons
        File? reference_amplicons_2 
        File? markersTable
        File? path_to_snv

        # Optional file for contamination detection pipeline
        File? barcodes_index

        File sample_metadata
        String out_prefix

        # TMP
        File panel_bed
        File CIGARVariants
        File ASVTable
        File ASVSeqs
        File ASV_to_CIGAR
        File ZeroReadsSampleList

	#REMOVE THIS AND PLACE IN SPECIFIC WDL        
        # Parameters for the contamination detection pipeline
        #Int minreads_threshold = 1000
        #Float contamination_threshold = 0.5
    }
    #FUTURE DEVELOPMENT
    # 1 CLEAN DOCKERFILES

#    call prepare_reference_files_t.prepare_reference_files as t_001_prepare_reference_files {
#        input:
#            #sample_ids = sample_ids,
#            panel_info = panel_info,
#            reference_genome = reference_genome,
#            reference_amplicons = reference_amplicons,
#            reference_amplicons_2 = reference_amplicons_2,
#            forward_primers_file = forward_primers_file,
#            reverse_primers_file = reverse_primers_file,
#            path_to_snv = path_to_snv
#    }
#
##    if(defined(barcodes_index)) {
##        call contamination_detection_t.contamination_detection {
##            input: 
##                #TMP FOR TESTING OF CONTAMINATION DETECTION PIPELINE
##                config_json = t_001_prepare_reference_files.config_json_out,
##                path_to_flist = t_001_prepare_reference_files.path_to_flist_o,
##                fastq1s = fastq1s,
##                fastq2s = fastq2s,
##                forward_primers_file = forward_primers_file,
##                reverse_primers_file = reverse_primers_file,
##                barcodes_index = barcodes_index
##        }
##    }
##
#
#    scatter (indx1 in range(length(fastq1s))) {
#        call cutadapters_t.cutadapters as t_002_cutadapters {
#            input:
#                fastq1 = fastq1s[indx1],
#                fastq2 = fastq2s[indx1]
#        }
#
#        call trimprimers_t.trimprimers as t_003_trimprimers {
#            input:
#                trimmed_fastq1 = t_002_cutadapters.fastq1_noadapters_o,
#                trimmed_fastq2 = t_002_cutadapters.fastq2_noadapters_o,
#                forward_primers_file = select_first([forward_primers_file, t_001_prepare_reference_files.forward_primers_o]),
#                reverse_primers_file = select_first([reverse_primers_file, t_001_prepare_reference_files.reverse_primers_o])
#        }
#    }
#
##    if(run_demultiplexing) {
##        call ampseq_pipeline_demult {
##            input: 
##                config_json = t_001_prepare_reference_files.config_json_out,
##                path_to_flist = path_to_flist,
##                fastq1s = fatsq1s,
##                fastq2s = fatsq2s,
##                forward_primers_file = forward_primers_file,
##                reverse_primers_file = reverse_primers_file,
##                reference = t_001_prepare_reference_files.reference_out,
##                reference2 = reference_amplicons_2,
##                path_to_snv = path_to_snv
##        }
##    }
#
#    call amplicon_denoising_t.amplicon_denoising as t_004_amplicon_denoising {
#        input:
#            sample_ids = sample_ids,
#            fastq1s = fastq1s,
#            fastq2s = fastq2s,
#            adaptor_rem1s = t_002_cutadapters.fastq1_noadapters_o,
#            adaptor_rem2s = t_002_cutadapters.fastq2_noadapters_o,
#            primer_rem1s = t_003_trimprimers.fastq1_no_primers_o,
#            primer_rem2s = t_003_trimprimers.fastq2_no_primers_o,
#            forward_primers_file = select_first([forward_primers_file, t_001_prepare_reference_files.forward_primers_o]),
#            reverse_primers_file = select_first([reverse_primers_file, t_001_prepare_reference_files.reverse_primers_o]),
#            reference_amplicons = t_001_prepare_reference_files.reference_o,
#            reference_amplicons_2 = reference_amplicons_2,
#            run_id = run_id,
#            path_to_snv = path_to_snv
#    }

    call asv_filtering_t.asv_filtering as t_005_asv_filtering{
        input:
            out_prefix = out_prefix,
            sample_metadata = sample_metadata,
            panel_bedfile = panel_bed,
            reference_amplicons = reference_amplicons,
            markersTable = markersTable,        
            reference_genome = reference_genome,
            CIGARVariants = CIGARVariants,
            ASVTable = ASVTable,
            ASVSeqs = ASVSeqs,
            ASV_to_CIGAR = ASV_to_CIGAR,
            ZeroReadsSampleList = ZeroReadsSampleList,
    }

    output {
        # PREPARE REFERENCE FILES
        #File path_to_flist_out = t_001_prepare_reference_files.path_to_flist_o
        #File reference_out = t_001_prepare_reference_files.reference_o
        #File panel_bedfile_out = t_001_prepare_reference_files.panel_bedfile_o
        #File? forward_primers_out = t_001_prepare_reference_files.forward_primers_o
        #File? reverse_primers_out = t_001_prepare_reference_files.reverse_primers_o
        
        # CUTADAPTERS AND TRIMPRIMERS
        #Array[File] cutadaptersout_f_out = t_002_cutadapters.fastq1_noadapters
        #Array[File] cutadaptersout_r_out = t_002_cutadapters.fastq2_noadapters
        #Array[File] trimprimersout_f_out = t_003_trimprimers.fastq1_no_primers
        #Array[File] trimprimersout_r_out = t_003_trimprimers.fastq2_no_primers

        # DADA2
        #File seqtab_out = t_004_amplicon_denoising.seqtab_o
        #File ASVBimeras_f = t_004_amplicon_denoising.ASVBimeras_o

        # PostProc_DADA2
        # File ASVTable_f = t_004_amplicon_denoising.ASVTable_o
        # File ASVSeqs_f = t_004_amplicon_denoising.ASVSeqs_o

        # ASV_to_CIGAR
        #File CIGARVariants_Bfilter_f = t_004_amplicon_denoising.CIGARVariants_Bfilter_o
        #File ASV_to_CIGAR_f = t_004_amplicon_denoising.ASV_to_CIGAR_o
        #File ZeroReadsSampleList_f = t_004_amplicon_denoising.ZeroReadsSampleList_o

        # ASV Filtering
        File ampseq_object_f = t_005_asv_filtering.ampseq_object_o
        File markersTable_f = t_005_asv_filtering.markersTable_o
        
        # CONTAMINATION DETECTION
        # File? contamination_detection_missing_files_f = contamination_detection.missing_files
        # File? contamination_detection_sample_cards_f = contamination_detection.contamination_detection_sample_cards
        # File? contamination_detection_report_f = contamination_detection.contamination_detection_report

    }
}
