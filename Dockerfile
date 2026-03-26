
FROM rocker/tidyverse:latest

RUN Rscript -e "install.packages(c('BiocManager'))" && \
    Rscript -e "BiocManager::install(c('DESeq2','apeglm','GSEABase','Category','GOstats'))"

