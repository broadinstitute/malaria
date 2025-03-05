version 1.0

task cutadapters {
    input {
        File fastq1
        File fastq2
        Int trim_galore_qvalue = 5
        Int trim_galore_length = 20
        Int downsample_fraction = 40000
    }
    # Get the basename of the fastq files without any fastq extensions:
    # CORRECT THIS: ALL BASENAMES GET THE R1 VALUE, EVEN THE REVERSE READS FILES
    String basename = basename(basename(basename(basename(fastq1, ".fq.gz"), ".fastq"), ".fq"), ".fastq.gz")
    command <<<
        export TMPDIR=tmp
        #set -euxo pipefail
        
        mkdir Results

        ################################################################### 
        # Trim adapters from the fastq files using TrimGalore             #
        ################################################################### 

        echo "Removing adapters"

        touch Results/~{basename}_val_1.fq.gz
        touch Results/~{basename}_val_2.fq.gz
        
        if [ "~{downsample_fraction}" -ne 0 ]; then
            echo "Downsampling to ~{downsample_fraction} reads."
            seqtk sample -s100 ~{fastq1} ~{downsample_fraction} > ~{basename}_R1_downsampled.fastq
            seqtk sample -s100 ~{fastq2} ~{downsample_fraction} > ~{basename}_R2_downsampled.fastq
            FASTQ1="~{basename}_R1_downsampled.fastq"
            FASTQ2="~{basename}_R2_downsampled.fastq"
        else
            echo "No downsampling requested."
            FASTQ1="~{fastq1}"
            FASTQ2="~{fastq2}"
        fi

        trim_galore --paired --gzip \
        --quality ~{trim_galore_qvalue} \
        --length ~{trim_galore_length} \
        --output_dir Results \
        --basename ~{basename} \
        -j 4 \
        "$FASTQ1" "$FASTQ2"

        echo "Done removing adapters"
    >>>

    output {
        # Generate a list of trimmed forward and reverse files
        File fastq1_noadapters_o = "Results/~{basename}_val_1.fq.gz"
        File fastq2_noadapters_o = "Results/~{basename}_val_2.fq.gz"
    }

    runtime {
        cpu: 4  # 4 cores is the sweet spot according to https://github.com/FelixKrueger/TrimGalore/blob/master/Docs/Trim_Galore_User_Guide.md
        memory: "16 GiB"
        disks: "local-disk 10 HDD"
        bootDiskSizeGb: 10
        preemptible: 3
        maxRetries: 1
        docker: 'jorgeamaya/cutadapters:v_0_0_2'
    }
}
