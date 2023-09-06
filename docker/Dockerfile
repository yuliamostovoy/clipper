FROM continuumio/miniconda3

MAINTAINER Yulia Mostovoy

# copy other resources
COPY ./environment.yml /
RUN wget https://raw.githubusercontent.com/yuliamostovoy/bigclipper/main/scripts/bigclipper_getclusters.py && \
    wget https://raw.githubusercontent.com/yuliamostovoy/bigclipper/main/scripts/bigclipper_processbam.py

# install conda packages
RUN conda env create -f /environment.yml && conda clean -a
RUN echo "conda activate bigclipper" >> ~/.bashrc
ENV PATH=/opt/conda/envs/bigclipper/bin/:/root/google-cloud-sdk/bin/:${PATH}

# install gsutil
RUN apt update
RUN apt install -y curl git-lfs parallel
RUN curl https://sdk.cloud.google.com | bash
