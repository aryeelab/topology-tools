FROM ubuntu:22.04
ARG GIT_TAG=__specify_in_build_arg__

RUN ln -fs /usr/share/zoneinfo/America/New_York /etc/localtime

RUN apt-get update && \
	apt-get -y install  bedtools \
						build-essential \
						git \
                        libbz2-dev  \
                        libcurl4-openssl-dev \
                        liblzma-dev \
                        libssl-dev \
						python3-deeptools \
						python3-markupsafe \	
						python3-pairtools \
						python3-pandas \
						python3-pillow \
						python3-pysam \
						python3-scipy \
						python3-tabulate \
						samtools \
						wget \
                        zlib1g-dev
												
# Install bwa
RUN	cd /tmp && \
	git clone https://github.com/lh3/bwa.git && \
	cd bwa; make && \
	cp bwa qualfa2fq.pl xa2multi.pl /usr/local/bin && \
	rm -fr /tmp/bwa

# Install preseq
RUN     cd /tmp && \
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
    
RUN echo $GIT_TAG > /VERSION