FROM condaforge/mambaforge:23.3.1-0

RUN mamba config \
    --add channels defaults \
    --add channels bioconda \
    --add channels conda-forge && \
    mamba create -n fileprep_env python=3.7.10 -y && \
    mamba install -n fileprep_env \
    bedtools=2.31.1 \
    pandas=1.3.0 \
    -c conda-forge -c bioconda && \
    mamba clean --all -f -y && \
    echo "source activate fileprep_env" > ~/.bashrc

ENV PATH /opt/conda/envs/fileprep_env/bin:$PATH
ENV PATH /opt/conda/envs/fileprep_env/bin/python:$PATH
SHELL ["conda", "run", "-n", "fileprep_env", "/bin/bash", "-c"]

COPY Code/add_entry_to_json.py Code/
COPY Code/check_input_headers.py Code/
