FROM debian:stable-20190228-slim
ENV DEBIAN_FRONTEND noninteractive
ENV MINICONDA_VERSION="4.5.11"
ENV CONDA_ENV_NAME="align_recalibrate"

ENV LANG=C.UTF-8 LC_ALL=C.UTF-8

#install packages like in continuumio/miniconda2:4.5.12
RUN apt-get update --fix-missing && apt-get install -y wget bzip2 ca-certificates \
    libglib2.0-0 libxext6 libsm6 libxrender1 \
    vim git mercurial subversion \
    gmsl gawk perl python \
    cython python-setuptools

#instal packages advised by Gianluca Della Vedova <gianluca@dellavedova.org>
RUN apt-get install -y procps uuid-runtime

WORKDIR /tmp
SHELL ["/bin/bash", "-c","-l"]

RUN wget -O miniconda.sh https://repo.anaconda.com/miniconda/Miniconda2-${MINICONDA_VERSION}-Linux-x86_64.sh 
RUN chmod ugo+x miniconda.sh 

RUN ./miniconda.sh -b -p /opt/conda && \
    rm miniconda.sh && \
    ln -s /opt/conda/etc/profile.d/conda.sh /etc/profile.d/conda.sh && \
    echo ". /opt/conda/etc/profile.d/conda.sh" >> /etc/profile


RUN conda config --add channels defaults && conda config --add channels bioconda && conda config --add channels conda-forge
RUN conda create -n ${CONDA_ENV_NAME} bcftools=1.9 bedtools=2.27 picard=2.18.15 samtools=1.9 trimmomatic=0.38  ucsc-liftover=357 bwa=0.7.17 fastqc=0.11.7 mosdepth=0.2.3 ensembl-vep=94.5 python=3.7.3

ENV PATH "/opt/conda/bin:/opt/conda/envs/${CONDA_ENV_NAME}/bin:${PATH}"
ENV CONDA_DEFAULT_ENV $CONDA_ENV_NAME
ENV CONDA_PREFIX /opt/conda/envs/$CONDA_ENV_NAME

# for interactive occam-run which uses -l and has the user home mounted by default
RUN echo "export PATH=/opt/conda/bin:/opt/conda/envs/${CONDA_ENV_NAME}/bin:${PATH}" >> /etc/profile

CMD ["/bin/bash"]
