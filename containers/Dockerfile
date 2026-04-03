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

RUN MAFFT_VERSION=7.520 && \
    wget -q https://mafft.cbrc.jp/alignment/software/mafft-${MAFFT_VERSION}-without-extensions-src.tgz && \
    tar -xzf mafft-${MAFFT_VERSION}-without-extensions-src.tgz && \
    cd mafft-${MAFFT_VERSION}-without-extensions/core && make install && \
    cd / && rm -rf mafft-${MAFFT_VERSION}-without-extensions*

RUN R --no-save --no-restore -e "\
    options(repos = c(CRAN = 'https://cloud.r-project.org')); \
    install.packages(c('tidyverse','data.table','here','vegan','changepoint', \
                       'magick','patchwork','gridExtra','ggrepel','fs', \
                       'scales','digest','tibble', \
                       'irr','ape','igraph','microeco'), dependencies = TRUE)"

RUN R --no-save --no-restore -e "\
    if (!requireNamespace('BiocManager', quietly = TRUE)) \
        install.packages('BiocManager', repos = 'https://cloud.r-project.org'); \
    BiocManager::install(c('Biostrings','DECIPHER'), ask = FALSE, update = FALSE)"

RUN R --no-save --no-restore -e "\
    pkgs <- c('tidyverse','Biostrings','DECIPHER','vegan','changepoint','magick','here', \
              'data.table','irr','ape','igraph','microeco'); \
    missing <- pkgs[!sapply(pkgs, requireNamespace, quietly = TRUE)]; \
    if (length(missing) > 0) stop(paste('Missing:', paste(missing, collapse=', '))); \
    message('All packages verified OK')"

RUN vsearch --version && mafft --version

LABEL org.opencontainers.image.title="SeCAT"
LABEL org.opencontainers.image.description="Sequence Consensus Amplicon Trimming v4.1"
LABEL org.opencontainers.image.version="4.1"
LABEL org.opencontainers.image.source="https://github.com/DerbyDT/SeCAT"
