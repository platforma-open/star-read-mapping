# Use R base image with specific version to match renv.lock files
FROM r-base:4.4.2

# Set environment variables
ENV R_VERSION=4.4.2
ENV BIOCONDUCTOR_VERSION=3.20
ENV RENV_VERSION=1.0.11

# Install system dependencies required for R packages
RUN apt-get update && apt-get install -y \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libcairo2-dev \
    libgit2-dev \
    default-libmysqlclient-dev \
    libpq-dev \
    libsasl2-dev \
    libsqlite3-dev \
    libssh2-1-dev \
    unixodbc-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    libgsl-dev \
    libbz2-dev \
    liblzma-dev \
    libpcre2-dev \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Install BiocManager and renv at system level
RUN R -e "install.packages(c('BiocManager', 'renv'), repos='https://packagemanager.posit.co/cran/latest')" \
    && R -e "BiocManager::install(version='$BIOCONDUCTOR_VERSION')"

# Set working directory
WORKDIR /app

# Copy project files (this will be specific to each R script when building)
COPY . ./

# Restore renv environment
RUN R -e "renv::consent(provided = TRUE)" && \
    R -e "renv::restore()"

# Default command runs Rscript
ENTRYPOINT ["Rscript"]
