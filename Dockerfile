FROM bioconductor/bioconductor_docker

RUN apt-get update && apt-get -y upgrade && \
apt-get install -y git && \
apt-get clean && apt-get purge && \
rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

RUN R -e "install.packages('optparse',dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('this.path',dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('ggplot2',dependencies=TRUE, repos='http://cran.rstudio.com/')"

RUN R -e "BiocManager::install('Rsamtools')"
RUN R -e "BiocManager::install('GenomicAlignments')"

ADD https://worldtimeapi.org/api/ip time.tmp
RUN git clone https://github.com/jlanej/long-read-plot.git
RUN Rscript long-read-plot/longReadPlot.R
CMD ["/bin/bash"]