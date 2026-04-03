# SeCAT v4.1
# Alignment: DECIPHER (R/Bioconductor)
# Clustering: VSEARCH
FROM rocker/r-ver:4.3.2

ENV LC_ALL=C
ENV PATH=/usr/local/bin:/usr/bin:/bin
ENV R_LIBS_USER=/usr/local/lib/R/library

# System dependencies
RUN apt-get update -qq && apt-get install -y --no-install-recommends \
        wget curl libcurl4-openssl-dev libssl-dev libxml2-dev \
        libfontconfig1-dev libharfbuzz-dev libfribidi-dev \
        libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev \
        libgit2-dev libmagick++-dev imagemagick build-essential \
        zlib1g-dev libbz2-dev liblzma-dev \
    && rm -rf /var/lib/apt/lists/*

# VSEARCH 2.23.0 — prebuilt binary
RUN VSEARCH_VERSION=2.23.0 && \
    wget -q https://github.com/torognes/vsearch/releases/download/v${VSEARCH_VERSION}/vsearch-${VSEARCH_VERSION}-linux-x86_64.tar.gz && \
    tar -xzf vsearch-${VSEARCH_VERSION}-linux-x86_64.tar.gz && \
    cp vsearch-${VSEARCH_VERSION}-linux-x86_64/bin/vsearch /usr/local/bin/vsearch && \
    chmod +x /usr/local/bin/vsearch && \
    rm -rf vsearch-${VSEARCH_VERSION}-linux-x86_64*

# Core CRAN packages
RUN R --no-save --no-restore -e "\
    options(repos = c(CRAN = 'https://cloud.r-project.org')); \
    install.packages(c('tidyverse','data.table','here','vegan','changepoint', \
                       'patchwork','gridExtra','ggrepel','fs', \
                       'scales','digest','tibble'), dependencies = TRUE)"

# magick — separate layer as it requires ImageMagick system libs
RUN R --no-save --no-restore -e "\
    options(repos = c(CRAN = 'https://cloud.r-project.org')); \
    install.packages('magick', dependencies = TRUE)"

# Ecological/network packages — separate layer as microeco has many deps
RUN R --no-save --no-restore -e "\
    options(repos = c(CRAN = 'https://cloud.r-project.org')); \
    install.packages(c('irr','ape','igraph'), dependencies = TRUE)"

RUN R --no-save --no-restore -e "\
    options(repos = c(CRAN = 'https://cloud.r-project.org')); \
    install.packages('microeco', dependencies = TRUE)"

# Bioconductor: Biostrings and DECIPHER
RUN R --no-save --no-restore -e "\
    if (!requireNamespace('BiocManager', quietly = TRUE)) \
        install.packages('BiocManager', repos = 'https://cloud.r-project.org'); \
    BiocManager::install(c('Biostrings','DECIPHER'), ask = FALSE, update = FALSE)"

# Verify — prints which packages are missing before failing
RUN R --no-save --no-restore -e "\
    pkgs <- c('tidyverse','Biostrings','DECIPHER','vegan','changepoint','magick','here', \
              'data.table','irr','ape','igraph','microeco'); \
    status <- sapply(pkgs, requireNamespace, quietly = TRUE); \
    cat('Package status:\n'); \
    for (i in seq_along(pkgs)) cat(sprintf('  %-20s %s\n', pkgs[i], ifelse(status[i], 'OK', 'MISSING'))); \
    missing <- pkgs[!status]; \
    if (length(missing) > 0) stop(paste('Missing packages:', paste(missing, collapse=', ')))"

RUN vsearch --version

LABEL org.opencontainers.image.title="SeCAT"
LABEL org.opencontainers.image.description="Sequence Consensus Amplicon Trimming v4.1"
LABEL org.opencontainers.image.version="4.1"
LABEL org.opencontainers.image.source="https://github.com/DerbyDT/SeCAT"
