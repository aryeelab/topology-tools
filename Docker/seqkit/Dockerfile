FROM ubuntu:22.04

RUN apt-get update && \
	apt-get -y install wget curl

RUN wget https://github.com/shenwei356/seqkit/releases/download/v2.3.1/seqkit_linux_amd64.tar.gz && \
	tar -zxvf *.tar.gz && \
	cp seqkit /usr/local/bin/
