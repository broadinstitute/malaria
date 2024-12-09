FROM condaforge/mambaforge:23.3.1-0

RUN mamba config \
    --add channels defaults \
    --add channels bioconda \
    --add channels conda-forge && \
    mamba create -n ampseq_env python=3.7.10 -y && \
    mamba install -n ampseq_env \
    r-base=4.1 \
    biopython=1.79 \
    muscle=3.8.1551 \
    bioconductor-dada2=1.20.0 \
    bioconductor-limma=3.48.0 \
    bioconductor-biostrings=2.60.0 \
    r-data.table=1.14.0 \
    r-viridisLite=0.4.0 \
    r-argparse=2.0.3 \
    r-seqinr=4.2_5 \
    r-stringdist=0.9.8 \
    -c conda-forge -c bioconda && \
    mamba clean --all -f -y && \
    echo "source activate ampseq_env" > ~/.bashrc

ENV PATH /opt/conda/envs/ampseq_env/bin:$PATH
ENV PATH /opt/conda/envs/ampseq_env/bin/python:$PATH
SHELL ["conda", "run", "-n", "ampseq_env", "/bin/bash", "-c"]

#    r-biocmanager=1.30.21 \
#    r-ape=5.7_1 \
#    r-vegan=2.6_4 \
#    r-rcolorbrewer=1.1_3 \
#    r-xlconnect=1.0.7 \
#    r-rmarkdown=2.22 \
#    r-gridextra=2.3 \
#    r-gt=0.9.0 \
#    r-tidyverse=2.0.0 \
#RUN apt-get update -y && apt-get install -y curl gnupg1 python3
#
#RUN echo "deb [signed-by=/usr/share/keyrings/cloud.google.gpg] http://packages.cloud.google.com/apt cloud-sdk main" | tee -a /etc/apt/sources.list.d/google-cloud-sdk.list && curl https://packages.cloud.google.com/apt/doc/apt-key.gpg | apt-key --keyring /usr/share/keyrings/cloud.google.gpg  add - && apt-get update -y && apt-get install google-cloud-sdk -y
#
#RUN curl https://sdk.cloud.google.com | bash
#ENV PATH=/root/google-cloud-sdk/bin/:${PATH}
#
#RUN apt-get install -y gcc python3-dev python3-setuptools && pip3 uninstall -y crcmod && pip3 install --no-cache-dir -U crcmod

COPY Code Code
