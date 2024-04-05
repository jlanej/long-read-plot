FROM ubuntu:22.04
# 22.04, jammy-20240227, jammy, latest
ARG DEBIAN_FRONTEND=noninteractive
ENV TZ=Etc/UTC

RUN apt-get update && apt-get -y upgrade && \
apt-get install -y wget

RUN wget -qO- https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc | sudo gpg --dearmor -o /usr/share/keyrings/r-project.gpg
RUN echo "deb [signed-by=/usr/share/keyrings/r-project.gpg] https://cloud.r-project.org/bin/linux/ubuntu jammy-cran40/" | sudo tee -a /etc/apt/sources.list.d/r-project.list

RUN apt-get update && apt-get -y upgrade && \
apt-get install -y libcurl4-openssl-dev git r-base r-base-dev

RUN R -e "sessionInfo();installed.packages()"

RUN apt-get update && apt-get -y upgrade && \
apt-get clean && apt-get purge && \
rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

RUN R -e "install.packages('optparse',dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('this.path',dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('ggplot2',dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('BiocManager',dependencies=TRUE, repos='http://cran.rstudio.com/')"

RUN R -e "BiocManager::install('Rsamtools')"
RUN R -e "BiocManager::install('GenomicAlignments')"

ADD https://worldtimeapi.org/api/ip time.tmp
RUN git clone https://github.com/jlanej/long-read-plot.git
RUN Rscript long-read-plot/longReadPlot.R
RUN R -e "sessionInfo();installed.packages()"
CMD ["/bin/bash"]