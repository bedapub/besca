FROM condaforge/mambaforge:4.12.0-0

LABEL MAINTAINER="paul.geser@roche.com"

SHELL ["/bin/bash", "-c"] 

RUN conda init bash
RUN /bin/bash -c "source ~/.bashrc"
RUN mkdir besca_base
RUN cd besca_base
COPY . . 
RUN mamba env create -f environment.lock.yml
RUN /bin/bash -c "source activate besca_dev"
