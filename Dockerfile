FROM ubuntu:22.04

RUN apt update && \
 apt install -y python3 python3-setuptools python3-pip git lz4

COPY . /genotype-error_src

RUN cd /genotype-error_src && pip3 install pyigd 

FROM ubuntu:22.04

COPY --from=0 /decompress.sh /usr/local/bin/decompress

RUN apt update && \
  apt install -y python3 python3-setuptools python3-pip git lz4 && \
  chmod +x /usr/local/bin/decompress
