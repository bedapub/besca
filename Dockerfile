FROM condaforge/mambaforge:4.12.0-0

LABEL MAINTAINER="paul.geser@roche.com"

SHELL ["/bin/bash", "-c"] 

RUN conda init bash
RUN source ~/.bashrc
RUN mkdir besca_base
RUN cd besca_base
COPY . . 
RUN mamba env create -f environment.lock.yml
RUN source activate besca_dev
