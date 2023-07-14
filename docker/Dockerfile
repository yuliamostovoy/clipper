FROM continuumio/miniconda3

MAINTAINER Yulia Mostovoy

# copy other resources
COPY ./environment.yml /
RUN wget https://github.com/yuliamostovoy/clipper/raw/master/split_reads.py && \
    wget https://github.com/yuliamostovoy/clipper/raw/master/split_reads_process_clusters.py

# install conda packages
RUN conda env create -f /environment.yml && conda clean -a
RUN echo "conda activate clipper" >> ~/.bashrc
ENV PATH=/opt/conda/envs/clipper/bin/:/root/google-cloud-sdk/bin/:${PATH}

# install gsutil
RUN apt update
RUN apt install -y curl git-lfs parallel
RUN curl https://sdk.cloud.google.com | bash
