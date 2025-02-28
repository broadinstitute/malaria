FROM condaforge/mambaforge:23.3.1-0

RUN mamba config \
    --add channels defaults \
    --add channels bioconda \
    --add channels conda-forge && \
    mamba create -n layout_env python=3.11.0 -y && \
    mamba install -n layout_env \
    r-base=4.2.3 \
    r-ade4=1.7_22 \
    r-dplyr=1.1.4 \
    r-magrittr=2.0.3 \
    r-tidyr=1.3.0 \
    r-ggplot2=3.4.4 \
    r-cowplot=1.1.2 \
    r-vegan=2.6_4 \
    r-ape=5.7_1 \
    r-rcolorbrewer=1.1_3 \
    r-hmisc=5.1_1 \
    r-ggpubr=0.6.0 \
    r-domc=1.3.8 \
    r-svmisc=1.2.3 \
    r-biocmanager=1.30.22 \
    bioconductor-biostrings=2.66.0 \
    r-dt=0.32 \
    r-stringr=1.5.0 \
    r-tmap=3.3_4 \
    r-sp=2.1_3 \
    r-argparse=2.2.2 \
    r-xlconnect=1.0.10 \
    r-igraph=2.0.2 \
    r-xlconnect=1.0.10 \
    unzip=6.0 \
    -c conda-forge -c bioconda && \
    mamba clean --all -f -y && \
    echo "source activate layout_env" > ~/.bashrc

ENV PATH /opt/conda/envs/layout_env/bin:$PATH
ENV PATH /opt/conda/envs/layout_env/bin/python:$PATH
ENV PATH /opt/conda/envs/layout_env/bin/unzip:$PATH
SHELL ["conda", "run", "-n", "layout_env", "/bin/bash", "-c"]

RUN apt-get update -y && apt-get install -y curl gnupg1 python3

RUN echo "deb [signed-by=/usr/share/keyrings/cloud.google.gpg] http://packages.cloud.google.com/apt cloud-sdk main" | tee -a /etc/apt/sources.list.d/google-cloud-sdk.list && curl https://packages.cloud.google.com/apt/doc/apt-key.gpg | apt-key --keyring /usr/share/keyrings/cloud.google.gpg  add - && apt-get update -y && apt-get install google-cloud-sdk -y

RUN curl https://sdk.cloud.google.com | bash
ENV PATH=/root/google-cloud-sdk/bin/:${PATH}

RUN apt-get install -y gcc python3-dev python3-setuptools && pip3 uninstall -y crcmod && pip3 install --no-cache-dir -U crcmod

Copy Code Code
