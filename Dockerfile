FROM r-base

RUN R -e "install.packages('optparse',dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('Rsamtools',dependencies=TRUE, repos='http://cran.rstudio.com/')"
ENV BASE=/app/

WORKDIR ${BASE}
COPY ./ ${BASE}
CMD ["Rscript", "longReadPlot.R"]