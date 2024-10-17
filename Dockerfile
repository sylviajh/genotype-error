FROM ubuntu:22.04

RUN apt update && \
 apt install -y python3 python3-setuptools python3-pip git lz4

COPY . /genotype-error_src

RUN cd /genotype-error_src && pip3 install pyigd 

