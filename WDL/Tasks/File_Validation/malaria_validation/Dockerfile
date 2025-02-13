
FROM condaforge/mambaforge:23.3.1-0

ENV DEBIAN_FRONTEND=noninteractive
RUN mamba config \
    --add channels defaults \
    --add channels bioconda \
    --add channels conda-forge && \
    mamba create -n validate_inputs python=3.7.10 -y && \
    mamba install -n validate_inputs \
    fastq_utils=0.25.2 \ 
    seqkit=2.9.0 \
    -c conda-forge -c bioconda && \
    mamba clean --all -f -y && \
    echo "source activate validate_inputs" > ~/.bashrc

ENV PATH=/opt/conda/envs/validate_inputs/bin:$PATH
ENV PATH=/opt/conda/envs/validate_inputs/bin/python:$PATH
SHELL ["conda", "run", "-n", "validate_inputs", "/bin/bash", "-c"]