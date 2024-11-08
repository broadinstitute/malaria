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

        mkdir Results

        ################################################################### 
        # Trim primers from the reads using cutadapt                      #
        ################################################################### 
        echo "Removing primers"

        echo "Primer files ~{forward_primers_file} and ~{reverse_primers_file}"
        echo "Removing primers from ~{trimmed_fastq1} and ~{trimmed_fastq2}"

        # Setting cores to 0 to use all available cores
        cutadapt -g file:~{forward_primers_file} \
            -G file:~{reverse_primers_file} \
            -o Results/~{basename}_mixed_op_1.fq.gz \
            -p Results/~{basename}_mixed_op_2.fq.gz \
            --pair-adapters --action=trim \
            --discard-untrimmed \
            -j 0 \
            "~{trimmed_fastq1}" "~{trimmed_fastq2}"

        echo "Primer removal done"

    >>>

    output {
        File fastq1_no_primers_o = "Results/~{basename}_mixed_op_1.fq.gz"
        File fastq2_no_primers_o = "Results/~{basename}_mixed_op_2.fq.gz"
    }
    runtime {
        cpu: 4
        memory: "5 GiB"
        disks: "local-disk 10 HDD"
        bootDiskSizeGb: 10
        preemptible: 3
        maxRetries: 1
        docker: 'jorgeamaya/trimprimers:v_1_0_0'
    }
}
