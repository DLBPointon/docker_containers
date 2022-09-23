FROM python:3.9
ARG TAG=0.001

USER root
################## MAINTAINER ######################
MAINTAINER Damon-Lee Pointon dp24@sanger.ac.uk
LABEL org.opencontainers.image.licenses="MIT"

################## INSTALLATION ######################
COPY --chmod=0755 scripts/filter_blast.py /usr/local/bin/
ENV PYTHONPATH=/usr/lib/python3.9/site-packages
RUN /usr/local/bin/python -m pip install --upgrade pip
RUN pip install numpy
RUN pip install pandas
