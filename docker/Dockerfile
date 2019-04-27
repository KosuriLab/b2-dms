# Start with rocker/tidyverse base image
FROM rocker/verse:3.5.2

# Install extra *nix utils
RUN apt-get update \
    && apt-get install -y \
    parallel \
    wget \
    pigz \
    vim \
    less

RUN echo 'will cite' | parallel --citation 1> /dev/null 2> /dev/null

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# R Deps

# x11, mesa, glu1 are so we can install paletteer
# cargo is for gifski which is for knitr
RUN apt-get update \
    && apt-get install -y \
    cargo \
    libx11-dev \
    mesa-common-dev \
    libglu1-mesa-dev

# Install extra R packages
RUN install2.r --error \
    --deps TRUE \
    paletteer \
    cowplot \
    gridExtra \
    ggrepel \
    ggbeeswarm \
    factoextra \
    dbscan \
    knitr \
    hexbin

# Install github packages
RUN installGithub.r \
    --deps TRUE \
    DavisVaughan/furrr \
    jlmelville/uwot \
    drsimonj/corrr

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Local Programs
WORKDIR /usr/local

# BBMAP
RUN apt-get update \
    && apt-get install -y \
    openjdk-8-jre

ENV BBMAP_DIR /usr/local/bbmap
RUN set -eu; \
    wget -O bbmap.tar.gz 'https://sourceforge.net/projects/bbmap/files/BBMap_37.36.tar.gz/download' \
    && tar xzf bbmap.tar.gz \
    && rm bbmap.tar.gz

ENV PATH ${BBMAP_DIR}:${PATH}

# freesasa
RUN apt-get update \
    && apt-get install -y \
    libjson-c-dev

RUN set -eu; \
    wget 'https://freesasa.github.io/freesasa-2.0.3.tar.gz' \
    && tar xzf freesasa-2.0.3.tar.gz \
    && cd freesasa-2.0.3 \
    && ./configure \
    && make -j8 \
    && make install \
    && cd .. \
    && rm freesasa-2.0.3.tar.gz

# parasail
RUN apt-get update \
    && apt-get install -y \
    autoconf \
    autotools-dev \
    libtool

RUN set -eu; \
    wget -O parasail-2.4.tar.gz 'https://github.com/jeffdaily/parasail/archive/v2.4.tar.gz' \
    && tar xzf parasail-2.4.tar.gz \
    && cd parasail-2.4 \
    && autoreconf -fi \
    && ./configure \
    && make -j8 \
    && ln -s /usr/local/parasail-2.4/apps/parasail_aligner /usr/local/bin/parasail_aligner \
    && cd .. \
    && rm parasail-2.4.tar.gz

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Python Env
# Install miniconda to /miniconda


ENV LANG=C.UTF-8 LC_ALL=C.UTF-8
ENV PATH /opt/conda/bin:$PATH

RUN apt-get update --fix-missing && \
    apt-get install -y \
    bzip2 \
    ca-certificates \
    curl \
    git \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-4.5.11-Linux-x86_64.sh -O ~/miniconda.sh && \
    /bin/bash ~/miniconda.sh -b -p /opt/conda && \
    rm ~/miniconda.sh && \
    /opt/conda/bin/conda clean -tipsy && \
    ln -s /opt/conda/etc/profile.d/conda.sh /etc/profile.d/conda.sh && \
    echo ". /opt/conda/etc/profile.d/conda.sh" >> ~/.bashrc && \
    echo "conda activate base" >> ~/.bashrc

WORKDIR /opt/conda
RUN conda update -y conda \
    && conda install -y \
    numpy \
    seaborn \
    pandas \
    scikit-learn \
    ipython \
    jupyter

RUN conda install -y -c conda-forge umap-learn

# barcode mapper and scripts all run on python2
RUN conda create -y -n py27 python=2.7 \
    && conda install -y -c conda-forge python-levenshtein

WORKDIR /home/rstudio/
