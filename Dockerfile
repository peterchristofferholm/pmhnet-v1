FROM nvcr.io/nvidia/pytorch:21.11-py3

## python packages ############################################################

RUN conda config --append channels bioconda \
    && conda install --yes \
        optuna=2.10.0 \
        snakemake=6.11.0 \
        pyyaml=5.4.1 \
        psycopg2=2.9.2 \
        numpy=1.21.4 \
        shap=0.40.0 \
        pandas=1.3.4 \
        umap-learn \
        numba

RUN pip3 install \
    eir-dl==0.1.24a0 \
    pycox==0.2.3

## install r and r-packages ###################################################
## inspired by code/dockerfiles from https://www.rocker-project.org/

ENV DEBIAN_FRONTEND noninteractive
ENV TZ=Europe/Copenhagen

ENV R_BASE_VERSION 4.2.1

RUN apt-get update -qq \
    && apt-get install --yes --no-install-recommends \
        software-properties-common dirmngr \  
    && wget -qO- https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc \
        | tee -a /etc/apt/trusted.gpg.d/cran_ubuntu_key.asc \
    && add-apt-repository "deb https://cloud.r-project.org/bin/linux/ubuntu $(lsb_release -cs)-cran40/" \
    && apt-get update -qq

RUN apt-get install -qq \
    r-base=${R_BASE_VERSION}-* \
    r-base-dev=${R_BASE_VERSION}-* \
    r-base-core=${R_BASE_VERSION}-* 

RUN apt-get install -qq \
        r-cran-littler littler \
    && ln -s /usr/lib/R/site-library/littler/examples/install.r /usr/local/bin/install.r \
    && ln -s /usr/lib/R/site-library/littler/examples/install2.r /usr/local/bin/install2.r \
    && ln -s /usr/lib/R/site-library/littler/examples/installBioc.r /usr/local/bin/installBioc.r \
    && ln -s /usr/lib/R/site-library/littler/examples/installDeps.r /usr/local/bin/installDeps.r \
    && ln -s /usr/lib/R/site-library/littler/examples/installGithub.r /usr/local/bin/installGithub.r \
    && ln -s /usr/lib/R/site-library/littler/examples/testInstalled.r /usr/local/bin/testInstalled.r \
    && install.r docopt

### tidyverse + tidymodels ###
RUN add-apt-repository ppa:c2d4u.team/c2d4u4.0+ \
    && apt-get update -qq \
    && apt-get install -qq \
        r-cran-tidyverse=1.3.1-* \
        r-cran-tidymodels=0.2.0-*

### database stuff ###
RUN apt-get install -qq \
    r-cran-odbc \
    r-cran-dbplyr \
    r-cran-rsqlite \
    unixodbc \
    unixodbc-dev \
    odbc-postgresql \
    libpq-dev \
    r-cran-rmariadb \
    r-cran-rms

COPY odbcinst.ini /etc/odbcinst.ini

RUN apt-get update -qq \
    && apt-get install -qq postgresql-client

### visualization stuff ###

RUN apt-get install -qq \
    r-cran-ggtext \
    r-cran-ggrepel \
    r-cran-ggforce \
    r-cran-patchwork \
    r-cran-ggdist \
    r-cran-ggridges

RUN add-apt-repository ppa:ubuntugis/ubuntugis-unstable \
    && apt-get update -qq  \
    && apt-get install --yes \
        libeigen3-dev \
        libv8-dev \
        libudunits2-dev \
        r-cran-sf \
        r-cran-rms \
        libjpeg-dev \
        libjpeg9

RUN install2.r --error --deps TRUE --skipinstalled --ncpus 8 \
    gggenes \
    geomtextpath \
    ggdensity \
    ggrastr \
    ggside 

RUN apt-get install -qq r-cran-remotes \
    && installGithub.r --deps TRUE \
        "hrbrmstr/ggalt" \
        "eclarke/ggbeeswarm"

### misc stuff ###

RUN apt-get install -qq \
    r-cran-yaml \
    r-cran-furrr \
    r-cran-tdroc \
    r-cran-vctrs \
    r-cran-rcppeigen \
    r-cran-zoo \
    r-cran-riskregression \ 
    r-cran-pec \ 
    r-cran-prodlim \ 
    r-cran-timeroc \ 
    r-cran-survivalroc \ 
    r-cran-timereg

RUN installGithub.r --deps TRUE "cboettig/arkdb"
