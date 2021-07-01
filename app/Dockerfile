FROM rocker/r-ver:3.6.3

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
    libjpeg62-turbo-dev \
    liblzma-dev \
    default-jre

# install conda & picrust2 (give shiny-user access to conda after installing shiny server further down)
RUN cd /opt/ && wget https://repo.anaconda.com/archive/Anaconda3-2020.11-Linux-x86_64.sh && bash /opt/Anaconda3-2020.11-Linux-x86_64.sh -b -p /opt/anaconda3
ENV PATH=/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/opt/anaconda3/bin
RUN conda create -n picrust2 -c bioconda -c conda-forge picrust2=2.3.0_b

# install cutadapt for primer trimming 
RUN conda create -n cutadapt -c bioconda cutadapt

# install fastqQC 
RUN cd /opt/ && wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.9.zip && unzip fastqc_v0.11.9.zip && chmod 777 /opt/FastQC/fastqc
ENV PATH=/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/opt/anaconda3/bin:/opt/FastQC

# Download and install shiny server
RUN wget --no-verbose https://download3.rstudio.org/ubuntu-14.04/x86_64/VERSION -O "version.txt" && \
    VERSION=$(cat version.txt)  && \
    wget --no-verbose "https://download3.rstudio.org/ubuntu-14.04/x86_64/shiny-server-$VERSION-amd64.deb" -O ss-latest.deb && \
    gdebi -n ss-latest.deb && \
    rm -f version.txt ss-latest.deb && \
    . /etc/environment && \
    R -e "install.packages(c('shiny', 'rmarkdown'), repos='$MRAN')" && \
    chown shiny:shiny /var/lib/shiny-server

ENV RENV_VERSION 0.13.1
RUN R -e "install.packages(c('remotes','huge'), repos = c(CRAN = 'https://cloud.r-project.org'))"
RUN R -e "install.packages('https://cran.r-project.org/src/contrib/Archive/foreign/foreign_0.8-76.tar.gz')"
RUN R -e "remotes::install_github('rstudio/renv@${RENV_VERSION}')"

COPY shiny-server.conf  /etc/shiny-server/shiny-server.conf
COPY shiny-server.sh /usr/bin/shiny-server.sh
COPY /R /srv/shiny-server
COPY ../renv.lock renv.lock


# download silva taxonomy reference
RUN wget https://zenodo.org/record/3986799/files/silva_nr99_v138_train_set.fa.gz && mv silva_nr99_v138_train_set.fa.gz /srv/shiny-server/data/taxonomy_annotation.fa.gz

RUN chown -R shiny:shiny /opt/anaconda3/*
RUN chown -R shiny:shiny /opt/FastQC/*
RUN chown -R shiny:shiny /srv/shiny-server
RUN R -e "renv::restore()"

EXPOSE 3838

CMD ["/usr/bin/shiny-server.sh"]

#sudo docker build -t gitlab.lrz.de:5005/biomedbigdata/namco:dev .
#sudo docker push gitlab.lrz.de:5005/biomedbigdata/namco:dev
#sudo docker run -d --name namco --rm -p 3838:3838 -v /srv/shinylog/:/var/log/shiny-server/ gitlab.lrz.de:5005/biomedbigdata/namco:latest

## docker image to gitlab:
#https://gitlab.lrz.de/help/user/packages/container_registry/index#build-and-push-images-from-your-local-machine
