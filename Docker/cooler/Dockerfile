FROM continuumio/miniconda

RUN apt-get update --fix-missing && \
    apt-get install -y unzip 
    
RUN conda install -y -c conda-forge -c bioconda cooler

COPY *.sizes /annotation/