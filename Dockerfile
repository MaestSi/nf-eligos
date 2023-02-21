FROM continuumio/miniconda3

MAINTAINER Simone Maestri <simone.maestri@iit.it>

RUN apt-get update  && \
    DEBIAN_FRONTEND=noninteractive \
    apt-get install -y \
    build-essential \
    wget \
    unzip \
    bzip2 \
    git \
    libidn11* \
    nano \
    less \
    wget \
    r-base \
    samtools \
    libcurl4-openssl-dev \
    libxml2-dev \
    libssl-dev \
    gfortran \
 && apt-get clean \
 && rm -rf /var/lib/apt/lists/*

RUN git clone https://gitlab.com/piroonj/eligos2.git 

# Create the environment:
RUN conda install -y -c bioconda -c conda-forge -c anaconda \
minimap2 python=3.6 pysam=0.13 pandas=0.23.4 pybedtools=0.8.0 bedtools=2.25 rpy2=2.8.5 r-base=3.4.1 tqdm=4.40.2 numpy=1.11.3

RUN R -e "install.packages('samplesizeCMH', repos='https://cloud.r-project.org')"

ENV PATH /eligos2:/eligos2/Scripts:$PATH
