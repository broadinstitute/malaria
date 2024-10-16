# Use an official Python image that is lightweight as the base
FROM python:3.8-slim

# Set environment variables to non-interactive mode to avoid prompts during install
ENV DEBIAN_FRONTEND=noninteractive

# Install dependencies and Perl
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
    perl \
    build-essential \
    gcc \
    libbz2-dev \
    zlib1g-dev \
    curl \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Install cutadapt using pip
RUN pip install --no-cache-dir cutadapt
