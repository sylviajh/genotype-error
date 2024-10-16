FROM ubuntu:22.04

RUN apt install -y git python3 lz4

COPY . /genotype-error_src

RUN cd /genotype-error_src && pip3 install pyigd 