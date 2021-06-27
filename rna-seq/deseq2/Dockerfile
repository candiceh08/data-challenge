FROM rticode/deseq2:1.26.0_08a6163
LABEL authors="candice.hermant@gmail.com" \
      description="Docker image containing all requirements for analysis of rna-seq gene-level count data"

# Install R packages

RUN R -e "install.packages('pheatmap', dependencies=TRUE, repos='https://cloud.r-project.org/')"
RUN R -e "install.packages('mclust', dependencies=TRUE, repos='https://cloud.r-project.org/')"
RUN R -e "install.packages('optparse', dependencies=TRUE, repos='https://cloud.r-project.org/')"
RUN R -e "install.packages('factoextra', dependencies=TRUE, repos='https://cloud.r-project.org/')"
RUN R -e "install.packages('cowplot', dependencies=TRUE, repos='https://cloud.r-project.org/')"
RUN R -e "install.packages('gridExtra', dependencies=TRUE, repos='https://cloud.r-project.org/')"
RUN R -e "install.packages('dplyr', dependencies=TRUE, repos='https://cloud.r-project.org/')"
RUN R -e "install.packages('stringr', dependencies=TRUE, repos='https://cloud.r-project.org/')"
RUN R -e "install.packages('openxlsx', dependencies=TRUE, repos='https://cloud.r-project.org/')"
RUN R -e "install.packages('RColorBrewer', dependencies=TRUE, repos='https://cloud.r-project.org/')"
RUN R -e "install.packages('colorspace', dependencies=TRUE, repos='https://cloud.r-project.org/')"
RUN R -e "install.packages('data.table', dependencies=TRUE, repos='https://cloud.r-project.org/')"