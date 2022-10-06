FROM python:latest

COPY qc_stats.py /home/

RUN pip install matplotlib
RUN pip install tabulate
RUN  apt-get update -y && \
     apt-get upgrade -y && \
     apt-get dist-upgrade -y && \
     apt-get -y autoremove && \
     apt-get clean
RUN apt-get install -y zip
