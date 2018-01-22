FROM debian:stretch
ENV LANG=C.UTF-8 LC_ALL=C.UTF-8

# Install procps (to get free and top)
RUN apt-get update --fix-missing && \
    apt-get install -y procps

# Install miniconda. Copied from https://hub.docker.com/r/continuumio/miniconda/~/dockerfile/ (because we need debian:stretch to install gcc 6 which is needed for HiCPro)
RUN apt-get install -y wget bzip2 ca-certificates \
    libglib2.0-0 libxext6 libsm6 libxrender1 \
    git mercurial subversion

RUN echo 'export PATH=/opt/conda/bin:$PATH' > /etc/profile.d/conda.sh && \
    wget --quiet https://repo.continuum.io/miniconda/Miniconda2-4.3.27-Linux-x86_64.sh -O ~/miniconda.sh && \
    /bin/bash ~/miniconda.sh -b -p /opt/conda && \
    rm ~/miniconda.sh

RUN apt-get install -y curl grep sed dpkg && \
    TINI_VERSION=`curl https://github.com/krallin/tini/releases/latest | grep -o "/v.*\"" | sed 's:^..\(.*\).$:\1:'` && \
    curl -L "https://github.com/krallin/tini/releases/download/v${TINI_VERSION}/tini_${TINI_VERSION}.deb" > tini.deb && \
    dpkg -i tini.deb && \
    rm tini.deb && \
    apt-get clean

ENV PATH /opt/conda/bin:$PATH

ENTRYPOINT [ "/usr/bin/tini", "--" ]
CMD [ "/bin/bash" ]
RUN apt-get update && apt-get install -y gcc && apt-get clean

RUN conda config --add channels bioconda
RUN conda install -c bioconda bx-python
RUN conda install -c bioconda scipy
RUN conda install -c bioconda numpy
RUN conda install -c bioconda pysam

# Install Bowtie2
RUN apt-get install -y gcc libtbb-dev unzip && \
        mkdir /src && \
        cd /src && \
        wget -O bowtie2.zip https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.3.2/bowtie2-2.3.2-linux-x86_64.zip/download && \
        unzip bowtie2.zip && \
        ln -s /src/bowtie2-2.3.2/bowtie2* /usr/local/bin


# Install R-3.4.0 (See https://cran.r-project.org/bin/linux/debian)
RUN apt-get install -y gnupg libssl-dev libcurl4-openssl-dev && \
        apt-key adv --keyserver keys.gnupg.net --recv-key 'E19F5F87128899B192B1A2C2AD5F960A256A04AF' && \
        echo "deb http://cloud.r-project.org/bin/linux/debian stretch-cran34/" >> /etc/apt/sources.list && \
        apt-get update && \
        apt-get install -y --force-yes r-base


#Installing necessary R packages
RUN Rscript -e "source('http://bioconductor.org/biocLite.R');biocLite('RColorBrewer');biocLite('ggplot2')"

# Install HiCPro
RUN cd /tmp && \
    wget -O HiC-Pro.tar.gz https://github.com/nservant/HiC-Pro/archive/v2.9.0.tar.gz && \
    tar zxf HiC-Pro.tar.gz && \
    cd HiC-Pro-2.9.0  && \ 
	sed 's/CLUSTER_SYS = /CLUSTER_SYS = LSF/g' config-install.txt  | sed 's/PREFIX =/PREFIX = \//g'  > tmp && \
	mv tmp config-install.txt && \
	ln -sf $(which gcc) /usr/local/bin/x86_64-conda_cos6-linux-gnu-gcc && \
	make configure && \
    make install && \
    cd .. && \
    rm -fr HiC-Pro* && \
    ln -s /HiC-Pro_2.9.0 /HiC-Pro

# Install a compatible version of gawk
RUN cd /tmp && \
    wget http://ftp.gnu.org/gnu/gawk/gawk-3.1.7.tar.gz && \
    tar zxf gawk-3.1.7.tar.gz && \
    cd gawk-3.1.7 && \
    ./configure --prefix=/usr --libexecdir=/usr/lib && \
    make && \
    make install && \
    cd .. && \
    rm -fr gawk-3.1.7* && \
    ln -sf /usr/bin/gawk /etc/alternatives/awk

# Install juicer tools
RUN apt-get install -y default-jre && \
    mkdir -p /usr/local/juicer && \
    cd /usr/local/juicer && \
    wget http://hicfiles.tc4ga.com.s3.amazonaws.com/public/juicer/juicer_tools.1.7.6_jcuda.0.8.jar

# Install bc for summing input fastq sizes
RUN apt-get install -y bc

COPY Mboi_resfrag_hg19.bed.gz /HiC-Pro/annotation
COPY chrom_hg19.sizes /HiC-Pro/annotation

COPY Mboi_resfrag_grch38.bed.gz /HiC-Pro/annotation
COPY chrom_grch38.sizes /HiC-Pro/annotation

COPY Mboi_resfrag_hg38.bed.gz /HiC-Pro/annotation
COPY chrom_hg38.sizes /HiC-Pro/annotation

COPY Mboi_resfrag_mm10.bed.gz /HiC-Pro/annotation
COPY chrom_mm10.sizes /HiC-Pro/annotation

RUN cd /HiC-Pro/annotation && \
    gunzip *.bed.gz