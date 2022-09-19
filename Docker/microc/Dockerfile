FROM python:latest


RUN apt-get update && \
      apt-get -y install sudo


copy ./installDep.sh .

RUN chmod +x /installDep.sh
RUN /installDep.sh