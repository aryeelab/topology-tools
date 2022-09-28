FROM ubuntu:22.04

RUN apt-get update && \
	apt-get -y install default-jre \
					   wget
												
# Install Juicer Tools
RUN     cd /usr/local/bin && \
        wget https://s3.amazonaws.com/hicfiles.tc4ga.com/public/juicer/juicer_tools_1.22.01.jar
    
ARG GIT_TAG=__specify_in_build_arg__    
RUN echo $GIT_TAG > /VERSION