FROM rocker/r-ver:4.2.0

RUN apt-get update && apt-get install -y \
    sudo \
    gdebi-core \
    pandoc \
    pandoc-citeproc \
    libcurl4-gnutls-dev \
    libcairo2-dev \
    libxt-dev \
    xtail \
    wget \
    libssl-dev \
    libxml2-dev \
    build-essential\
    libv8-dev \
    vim \
    bowtie2 \
    bzip2 \
    libbz2-dev \
    libjpeg-dev \
    liblzma-dev \
    default-jre \
    libgsl-dev \
    cmake
    
    
# Download and install shiny server
RUN wget --no-verbose https://download3.rstudio.org/ubuntu-14.04/x86_64/VERSION -O "version.txt" && \
    VERSION=$(cat version.txt)  && \
    wget --no-verbose "https://download3.rstudio.org/ubuntu-14.04/x86_64/shiny-server-$VERSION-amd64.deb" -O ss-latest.deb && \
    gdebi -n ss-latest.deb && \
    rm -f version.txt ss-latest.deb && \
    . /etc/environment && \
    R -e "install.packages(c('shiny', 'rmarkdown'), repos='$MRAN')" && \
    chown shiny:shiny /var/lib/shiny-server


# install conda (give shiny-user access to conda after installing shiny server further down)
RUN cd /opt/ && wget https://repo.anaconda.com/miniconda/Miniconda3-py37_4.11.0-Linux-x86_64.sh && bash /opt/Miniconda3-py37_4.11.0-Linux-x86_64.sh -b -p /opt/miniconda3
ENV PATH=/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/opt/miniconda3/bin

# change permissions for conda
RUN chown -R shiny:shiny /opt/miniconda3/*

# new conda env with picrust2 and lotus2
COPY conda_env.yml /etc/conda_env.yml
RUN conda env create -f /etc/conda_env.yml
# download usearch binary
RUN cd /opt/ && wget http://www.drive5.com/downloads/usearch11.0.667_i86linux32.gz && gunzip usearch* && chmod 777 usearch* && mv usearch* usearch 
# link usearch binary to lotus2
RUN conda run -n namco_env lotus2 -link_usearch /opt/usearch

# clean up temporary conda files
RUN conda clean --all

# install fastqQC 
RUN cd /opt/ && wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.9.zip && unzip fastqc_v0.11.9.zip && chmod 777 /opt/FastQC/fastqc
ENV PATH=/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/opt/miniconda3/bin:/opt/FastQC

# download FastTree binary and set number of cores
RUN cd /opt/ && wget http://www.microbesonline.org/fasttree/FastTreeMP && chmod 777 /opt/FastTreeMP
ENV OMP_NUM_THREADS=3

# setup renv to handle R packages
ENV RENV_VERSION 0.15.2
RUN R -e "install.packages(c('remotes','huge'), repos = c(CRAN = 'https://cloud.r-project.org'))"
RUN R -e "install.packages('https://cran.r-project.org/src/contrib/Archive/foreign/foreign_0.8-76.tar.gz')"
RUN R -e "remotes::install_github('rstudio/renv@${RENV_VERSION}')"

# install all packages
COPY renv.lock renv.lock
RUN R -e "renv::restore()"

# copy code and configs into container
COPY shiny-server.conf  /etc/shiny-server/shiny-server.conf
COPY shiny-server.sh /usr/bin/shiny-server.sh
COPY /R /srv/shiny-server
RUN chown -R shiny:shiny /srv/shiny-server

# download silva taxonomy reference
RUN wget https://zenodo.org/record/3986799/files/silva_nr99_v138_train_set.fa.gz && mv silva_nr99_v138_train_set.fa.gz /srv/shiny-server/data/taxonomy_annotation.fa.gz


EXPOSE 3838

#CMD ["/usr/bin/shiny-server.sh"]
CMD ["R", "-e", "shiny::runApp('/srv/shiny-server/', host = '0.0.0.0', port = 3838)"]
