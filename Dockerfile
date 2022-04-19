FROM python:3.9.12-alpine

LABEL MAINTAINER="paul.geser@roche.com"


RUN pip install pytest urllib3 pandas=1.4.1 pkg_resources scanpy
RUN mkdir besca_base
RUN cd besca_base
COPY . . 
RUN cd besca/datasets
RUN python -c 'from _datasets import *; print Baron2016_raw()'

