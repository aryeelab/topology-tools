FROM ubuntu:22.04

RUN ln -fs /usr/share/zoneinfo/America/New_York /etc/localtime

RUN apt-get update && \
    apt-get -y install  bedtools \
                        tabix \
                        build-essential \
                        curl \
                        git \
                        libbz2-dev  \
                        libcurl4-openssl-dev \
                        liblzma-dev \
                        libssl-dev \
                        python3-pip \
                        samtools \
                        wget \
                        zlib1g-dev

# Install bwa
RUN cd /tmp && \
    git clone https://github.com/lh3/bwa.git && \
    cd bwa; make && \
    cp bwa qualfa2fq.pl xa2multi.pl /usr/local/bin && \
    rm -fr /tmp/bwa

# Install preseq
RUN cd /tmp && \
    wget https://github.com/samtools/htslib/releases/download/1.16/htslib-1.16.tar.bz2 && \
    tar jxf htslib-1.16.tar.bz2 && \
    cd htslib-1.16 && \
    ./configure && make && make install && \
    cd .. && \
    rm -fr htslib-1.16 && \
    wget https://github.com/smithlabcode/preseq/releases/download/v3.1.2/preseq-3.1.2.tar.gz && \
    tar zxf preseq-3.1.2.tar.gz && \
    cd preseq-3.1.2 && \
    ./configure --enable-hts && make && make install && \
    cd .. && \
    rm -fr preseq-3.1.2

# Install Python packages via pip to avoid numpy/scipy version conflicts
# Use numpy<1.24 to avoid np.int deprecation issues with pairtools
RUN pip3 install "numpy<1.24" scipy pandas pysam==0.19.1 deeptools

# Install pairtools
# We use a patched version of pairtools since Mapping has moved from collections to collections.abc
# (See, for example: https://github.com/tensorflow/tensorboard/issues/5478)
# Need to install cython first and ensure proper compilation of extensions
RUN pip3 install cython click nose
RUN cd /tmp && \
    git clone https://github.com/aryeelab/pairtools.git && \
    cd pairtools && \
    git checkout v0.3.0p && \
    python3 setup.py build_ext --inplace && \
    pip3 install .



# Install Google Cloud SDK (primarily for gsutil)
RUN curl -sSL https://sdk.cloud.google.com  | bash && \
    ln -s /root/google-cloud-sdk/bin/gcloud /usr/local/bin/ && \
    ln -s /root/google-cloud-sdk/bin/gsutil /usr/local/bin/

ARG GIT_TAG=__specify_in_build_arg__
RUN echo $GIT_TAG > /VERSION