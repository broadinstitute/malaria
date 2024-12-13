# Use an official Ubuntu image as a base
FROM ubuntu:20.04

# Set non-interactive mode to prevent issues with package installation
ENV DEBIAN_FRONTEND=noninteractive

# Install dependencies
RUN apt-get update && apt-get install -y \
    wget \
    unzip \
    curl \
    perl \
    python3 \
    python3-pip \
    git \
    ca-certificates \
    gzip \
    pigz \
    bzip2 \
    seqtk \
    && rm -rf /var/lib/apt/lists/*

# Install Cutadapt (required by Trim Galore)
RUN pip3 install cutadapt

# Install Trim Galore
RUN wget https://github.com/FelixKrueger/TrimGalore/archive/refs/tags/0.6.6.zip -O trim_galore.zip && \
    unzip trim_galore.zip && \
    mv TrimGalore-0.6.6 /opt/trim_galore && \
    ln -s /opt/trim_galore/trim_galore /usr/local/bin/trim_galore
