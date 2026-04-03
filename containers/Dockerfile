# SeCAT v4.1
FROM rocker/r-ver:4.3.2

ENV LC_ALL=C
ENV PATH=/usr/local/bin:/usr/bin:/bin
ENV R_LIBS_USER=/usr/local/lib/R/library

RUN apt-get update -qq && apt-get install -y --no-install-recommends \
        wget curl libcurl4-openssl-dev libssl-dev libxml2-dev \
        libfontconfig1-dev libharfbuzz-dev libfribidi-dev \
        libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev \
        libgit2-dev libmagick++-dev imagemagick build-essential \
        zlib1g-dev libbz2-dev liblzma-dev \
    && rm -rf /var/lib/apt/lists/*

RUN VSEARCH_VERSION=2.23.0 && \
    wget -q https://github.com/torognes/vsearch/releases/download/v${VSEARCH_VERSION}/vsearch-${VSEARCH_VERSION}-linux-x86_64.tar.gz && \
    tar -xzf vsearch-${VSEARCH_VERSION}-linux-x86_64.tar.gz && \
    cp vsearch-${VSEARCH_VERSION}-linux-x86_64/bin/vsearch /usr/local/bin/vsearch && \
    chmod +x /usr/local/bin/vsearch && \
    rm -rf vsearch-${VSEARCH_VERSION}-linux-x86_64*

RUN R --no-save --no-restore -e "\
    options(repos = c(CRAN = 'https://cloud.r-project.org')); \
    install.packages(c('ggplot2','dplyr','tidyr','readr','purrr','tibble', \
                       'stringr','forcats','lubridate','tidyverse'), \
                     dependencies = TRUE)"

RUN R --no-save --no-restore -e "\
    options(repos = c(CRAN = 'https://cloud.r-project.org')); \
    install.packages(c('data.table','here','vegan','changepoint', \
                       'magick','patchwork','gridExtra','ggrepel', \
                       'fs','scales','digest'), dependencies = TRUE)"

RUN R --no-save --no-restore -e "\
    if (!requireNamespace('BiocManager', quietly = TRUE)) \
        install.packages('BiocManager', repos = 'https://cloud.r-project.org'); \
    BiocManager::install(c('Biostrings','DECIPHER'), ask = FALSE, update = FALSE)"

RUN R --no-save --no-restore -e "\
    pkgs <- c('ggplot2','dplyr','tidyr','readr','Biostrings','DECIPHER', \
              'vegan','changepoint','magick','here','data.table'); \
    missing <- pkgs[!sapply(pkgs, requireNamespace, quietly = TRUE)]; \
    if (length(missing) > 0) stop(paste('Missing:', paste(missing, collapse=', '))); \
    message('All packages verified OK')"

RUN vsearch --version

LABEL org.opencontainers.image.title="SeCAT"
LABEL org.opencontainers.image.description="Sequence Consensus Amplicon Trimming v4.1"
LABEL org.opencontainers.image.version="4.1"
LABEL org.opencontainers.image.source="https://github.com/DerbyDT/SeCAT"
