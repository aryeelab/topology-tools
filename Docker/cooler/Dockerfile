FROM continuumio/miniconda
ENV LANG=C.UTF-8 LC_ALL=C.UTF-8

RUN conda install -y -c conda-forge -c bioconda cooler

COPY *.sizes /annotation/