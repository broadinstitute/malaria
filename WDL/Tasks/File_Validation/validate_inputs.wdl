version 1.0

task validate_fastqs {
    input {
        File fastq1
        File fastq2
        Boolean? use_original_files = true
    }

    String fq1_prefix = basename(sub(fastq1, "\.fq", "\.fastq"), ".fastq.gz")
    String fq2_prefix = basename(sub(fastq2, "\.fq", "\.fastq"), ".fastq.gz")
                                    
    command <<< 
    set -euo pipefail
    # Validate FASTQs
    fastq_info ~{fastq1} ~{fastq2}
    # Get basic statistics
    seqkit stats ~{fastq1} ~{fastq2}

    # Properly pair FASTQs if indicated
    if [[ "~{use_original_files}" == "true" ]]; then
        echo "Using original FASTQ files."
        cp ~{fastq1} ~{fq1_prefix}.fastq.gz
        cp ~{fastq2} ~{fq2_prefix}.fastq.gz
    else
        echo "Pairing FASTQ files."
        mkdir Paired
        seqkit pair -1 ~{fastq1} -2 ~{fastq2} -O Paired
        mv Paired/~{fq1_prefix}.fastq.gz ~{fq1_prefix}.fastq.gz
        mv Paired/~{fq2_prefix}.fastq.gz ~{fq2_prefix}.fastq.gz
    fi

    >>>
    output {
        File fastq1_o = fq1_prefix + ".fastq.gz"
        File fastq2_o = fq2_prefix + ".fastq.gz"
    }
    
    runtime {
        cpu: 1
        memory: "1 GB"
		disks: "local-disk 10 HDD"
		bootDiskSizeGb: 10
		preemptible: 3
		maxRetries: 1
		docker: "jorgeamaya/validate_inputs:v_1_0_0"
    }
}