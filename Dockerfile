FROM condaforge/mambaforge:4.12.0-0

LABEL MAINTAINER="paul.geser@roche.com"

RUN conda init bash
RUN source ~/.bashrc
RUN mamba env create -f environment.lock.yml
RUN source activate besca_dev
