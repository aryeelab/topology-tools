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
                        python3-markupsafe \    
                        python3-pandas \
                        python3-pillow \
                        python3-pip \
                        python3-scipy \
                        python3-tabulate \
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

# Install a recent versions of pysam 
# As of 2022-10-04 the apt repo version of pysam (0.17.0+ds-2build1) is too old for pairtools 1.0.1
RUN pip3 install pysam==0.19.1
    
# Install deeptools
RUN pip3 install deeptools 

# Install pairtools 
# We use a patched version of pairtools since Mapping has moved from collections to collections.abc
# (See, for example: https://github.com/tensorflow/tensorboard/issues/5478)
RUN cd /tmp && \
    git clone https://github.com/aryeelab/pairtools.git && \
    cd pairtools && \
    git checkout v0.3.0p && \
    pip3 install . && \
    pip3 install .

# Install pbgzip
RUN cd /tmp && \
	git clone https://github.com/nh13/pbgzip.git && \
	cd pbgzip && \
	git checkout 2b09f97b5f20b6d83c63a5c6b408d152e3982974 && \
	sh autogen.sh && \
	./configure && \
	make && \
	make install

# Install Google Cloud SDK (primarily for gsutil)
RUN curl -sSL https://sdk.cloud.google.com  | bash && \
    ln -s /root/google-cloud-sdk/bin/gcloud /usr/local/bin/ && \
    ln -s /root/google-cloud-sdk/bin/gsutil /usr/local/bin/

ARG GIT_TAG=__specify_in_build_arg__    
RUN echo $GIT_TAG > /VERSION