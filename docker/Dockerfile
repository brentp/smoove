FROM ubuntu:16.04

ENV PATH /opt/conda/bin:$PATH

RUN apt-get update --fix-missing && \
    apt-get install -qy wget curl git bzip2 ca-certificates procps zlib1g-dev \
    make build-essential cmake libncurses-dev ncurses-dev g++ gcc \
    nfs-common pigz bedtools gawk fuse mdadm time \
    libbz2-dev lzma-dev liblzma-dev libglib2.0-0 libxext6 libsm6 libxrender1 \
    syslog-ng libssl-dev libtool autoconf automake \
    libcurl4-openssl-dev libffi-dev libblas-dev liblapack-dev \
    libatlas-base-dev libroot-math-mathmore-dev

RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda2-4.5.11-Linux-x86_64.sh -O ~/miniconda.sh && \
    /bin/bash ~/miniconda.sh -b -p /opt/conda && \
    rm ~/miniconda.sh && \
    ln -s /opt/conda/etc/profile.d/conda.sh /etc/profile.d/conda.sh && \
    echo ". /opt/conda/etc/profile.d/conda.sh" >> ~/.bashrc && \
    echo "conda activate base" >> ~/.bashrc

RUN conda create -n smoove-env -c conda-forge -c bioconda python=2.7 awscli numpy scipy cython pysam toolshed pyvcf pyfaidx cyvcf2 svtyper svtools
RUN echo "source activate smoove-env" > ~/.bashrc
ENV PATH /opt/conda/envs/smoove-env/bin:$PATH
ENV HTSLIB_LIBRARY_DIR /usr/local/lib
ENV HTSLIB_INCLUDE_DIR /usr/local/include
ENV LD_LIBRARY_PATH /usr/local/lib

COPY ./docker-build.sh .
RUN bash docker-build.sh
WORKDIR /work/
