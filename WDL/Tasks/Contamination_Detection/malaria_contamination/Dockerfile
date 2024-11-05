FROM condaforge/mambaforge:23.3.1-0

RUN mamba config \
    --add channels defaults \
    --add channels bioconda \
    --add channels conda-forge && \
    mamba create -n ci_env python=3.7.10 -y && \
    mamba install -n ci_env \
    r-base=4.1 \
    bbmap=39.01 \
    biopython=1.79 \
    trim-galore=0.6.6 \
    cutadapt=3.4 \
    muscle=3.8.1551 \
    r-viridisLite=0.4.0 \
    r-argparse=2.0.3 \
    r-stringdist=0.9.8 \
    r-rmarkdown=2.22 \
    r-gridextra=2.3 \
    r-gt=0.9.0 \
    r-tidyverse=2.0.0 \
    r-reshape2=1.4.4 \
    r-plyr=1.8.6 \
    -c conda-forge -c bioconda && \
    mamba clean --all -f -y && \
    echo "source activate ci_env" > ~/.bashrc

ENV PATH /opt/conda/envs/ci_env/bin:$PATH
ENV PATH /opt/conda/envs/ci_env/bin/python:$PATH
SHELL ["conda", "run", "-n", "ci_env", "/bin/bash", "-c"]

COPY Code Code
