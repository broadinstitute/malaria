version 1.0

task cutadapters {
    input {
        File fastq1
        File fastq2
        Int trim_galore_qvalue
        Int trim_galore_length
    }

    # Get the basename of the fastq files without any fastq extensions:
    String basename = basename(basename(basename(fastq1, ".fq.gz"), ".fastq"), ".fq")

    command <<<
        export TMPDIR=tmp
        set -euxo pipefail
        
        echo "Removing adapters"
        mkdir Results
        
        trim_galore --paired --gzip \
            --quality ~{trim_galore_qvalue} \
            --length ~{trim_galore_length} \
            --output_dir Results \
            --basename ~{basename} \
            -j 4 \ 
            ~{fastq1} ~{fastq2}
    >>>

    output {
        # Generate a list of trimmed forward and reverse files
        File fastq1_noadapters = "Results/~{basename}_val_1.fq.gz"
        File fastq2_noadapters = "Results/~{basename}_val_2.fq.gz"
    }

    runtime {
        cpu: 4  # 4 cores is the sweet spot according to https://github.com/FelixKrueger/TrimGalore/blob/master/Docs/Trim_Galore_User_Guide.md
        memory: "16 GiB"
        disks: "local-disk 10 HDD"
        bootDiskSizeGb: 10
        preemptible: 3
        maxRetries: 1
        docker: 'jorgeamaya/cutadapters:testing'
    }
}
