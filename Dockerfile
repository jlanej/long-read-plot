FROM ubuntu:latest
USER root
ARG DEBIAN_FRONTEND=noninteractive
ENV TZ=Etc/UTC
RUN apt-get update && apt-get -y upgrade && \
apt-get install -y build-essential wget libncurses5-dev zlib1g-dev libbz2-dev liblzma-dev libcurl3-dev git  r-base r-base-dev && \
apt-get clean && apt-get purge && \
rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# 
# create a .Rprofile file in the home directory
# options(HTTPUserAgent = sprintf("R/%s R (%s)", getRversion(), paste(getRversion(), R.version["platform"], R.version["arch"], R.version["os"])))
RUN echo "options(HTTPUserAgent = sprintf(\"R/%s R (%s)\", getRversion(), paste(getRversion(), R.version[\"platform\"], R.version[\"arch\"], R.version[\"os\"])))" > $HOME/.Rprofile



RUN R -e "install.packages('crayon')"
RUN R -e "library('crayon')"
RUN R -e "source('https://docs.rstudio.com/rspm/admin/check-user-agent.R')"

RUN R -e "install.packages('optparse',dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('BiocManager',dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "BiocManager::install('Biostrings')"
RUN R -e "BiocManager::install('Rsamtools')"

ADD https://worldtimeapi.org/api/ip time.tmp

WORKDIR /app/
RUN git clone https://github.com/jlanej/long-read-plot.git
WORKDIR /app/long-read-plot/
RUN Rscript longReadPlot.R
CMD ["/bin/bash"]
