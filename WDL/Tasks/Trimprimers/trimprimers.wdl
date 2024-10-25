version 1.0

task trimprimers {
    input {
        File trimmed_fastq1
        File trimmed_fastq2
        File forward_primers_file
        File reverse_primers_file
    }

    String basename = basename(trimmed_fastq1, "_val_1.fq.gz")

    command <<<
        export TMPDIR=tmp
        set -euxo pipefail

        echo "Removing primers"
        mkdir Results

        cutadapt -g file:~{forward_primers_file} \
            -G file:~{reverse_primers_file} \
            -o Results/${sample}_mixed_op_1.fq.gz \
            -p Results/${sample}_mixed_op_2.fq.gz \
            --pair-adapters --action=trim \
            --discard-untrimmed \
            -j 0 \  # setting cores to 0 to use all available cores
            ~{trimmed_fastq1} ~{trimmed_fastq2}

    >>>
    output {
        File fastq1_no_primers = "Results/~{basename}_mixed_op_1.fq.gz"
        File fastq2_no_primers = "Results/~{basename}_mixed_op_2.fq.gz"
    }
    runtime {
        cpu: 4
        memory: "5 GiB"
        disks: "local-disk 10 HDD"
        bootDiskSizeGb: 10
        preemptible: 0
        maxRetries: 1
        docker: 'jorgeamaya/trimprimers:testing'
    }
}
