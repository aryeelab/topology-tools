FROM ubuntu:18.04

### Set conda environment
ENV PATH="/root/miniconda3/bin:${PATH}"
ARG PATH="/root/miniconda3/bin:${PATH}"

## install sudo and wget
RUN apt-get update && \
      apt-get -y install sudo

RUN apt-get install -y wget && rm -rf /var/lib/apt/lists/*

## Install anaconda
RUN wget \
    https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh \
    && mkdir /root/.conda \
    && bash Miniconda3-latest-Linux-x86_64.sh -b \
    && rm -f Miniconda3-latest-Linux-x86_64.sh 
RUN conda --version

## copy the script to install micro-c dependencies
copy ./installDep.sh .
RUN chmod +x /installDep.sh

## install system requirements
RUN sudo apt-get -y update && sudo apt-get -y upgrade

RUN apt-get install -y build-essential python3 python3-pip

RUN sudo update-alternatives --install /usr/bin/python python /usr/bin/python3 1
RUN sudo update-alternatives --install /usr/bin/pip pip /usr/bin/pip3 1

RUN sudo apt-get -y install bedtools

RUN sudo apt install -y git

## Install dependencies
RUN /installDep.sh