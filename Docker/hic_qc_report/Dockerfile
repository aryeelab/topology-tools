FROM debian:stretch
ENV LANG=C.UTF-8 LC_ALL=C.UTF-8

RUN apt-get update --fix-missing && \
    apt-get install -y unzip git

# Install R-3.4.0 (See https://cran.r-project.org/bin/linux/debian)
RUN apt-get install -y gnupg libssl-dev libcurl4-openssl-dev && \
        apt-key adv --keyserver keys.gnupg.net --recv-key 'E19F5F87128899B192B1A2C2AD5F960A256A04AF' && \
        echo "deb http://cloud.r-project.org/bin/linux/debian stretch-cran34/" >> /etc/apt/sources.list && \
        apt-get update && \
        apt-get install -y --force-yes r-base


# Installing necessary R packages
RUN Rscript -e "source('http://bioconductor.org/biocLite.R'); biocLite(c('RColorBrewer', 'RColorBrewer', 'ggplot2', 'dplyr', 'reshape2', 'cowplot', 'optparse'))"

COPY HiCStats.R /usr/local/bin/

