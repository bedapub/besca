FROM condaforge/mambaforge:4.12.0-0

LABEL MAINTAINER="paul.geser@roche.com"

RUN conda init bash
#RUN /bin/bash -c 'source /opt/ros/melodic/setup.bash'
RUN mamba env create -f environment.lock.yml
RUN source activate besca_dev
