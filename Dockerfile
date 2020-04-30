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
    build-essential
   

# Download and install shiny server
RUN wget --no-verbose https://download3.rstudio.org/ubuntu-14.04/x86_64/VERSION -O "version.txt" && \
    VERSION=$(cat version.txt)  && \
    wget --no-verbose "https://download3.rstudio.org/ubuntu-14.04/x86_64/shiny-server-$VERSION-amd64.deb" -O ss-latest.deb && \
    gdebi -n ss-latest.deb && \
    rm -f version.txt ss-latest.deb && \
    . /etc/environment && \
    R -e "install.packages(c('shiny', 'rmarkdown'), repos='$MRAN')" && \
    chown shiny:shiny /var/lib/shiny-server

ENV RENV_VERSION 0.9.3-69
RUN R -e "install.packages(c('remotes','huge'), repos = c(CRAN = 'https://cloud.r-project.org'))"
RUN R -e "remotes::install_github('rstudio/renv@${RENV_VERSION}')"

COPY shiny-server.conf  /etc/shiny-server/shiny-server.conf
COPY shiny-server.sh /usr/bin/shiny-server.sh
COPY /R /srv/shiny-server
COPY renv.lock renv.lock

RUN R -e "renv::restore()"

EXPOSE 3838

CMD ["/usr/bin/shiny-server.sh"]

#sudo docker build -t namco-app .
#sudo docker run -d --name namco --rm -p 3838:3838 -v /srv/shinylog/:/var/log/shiny-server/ namco-app 