FROM rocker/rstudio:latest

# Get pkgs
RUN apt-get update && apt-get install -y --no-install-recommends \
    pandoc git make build-essential wget ca-certificates gdebi-core \
    libcurl4-openssl-dev libssl-dev libxml2-dev libgit2-dev \
    libharfbuzz-dev libfribidi-dev libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev zlib1g-dev \
 && rm -rf /var/lib/apt/lists/*

# Core R packages you need
RUN R -q -e "install.packages(c('rmarkdown','knitr','targets','tarchetypes','tidyverse','data.table','janitor','here','glue','whisker','devtools'), repos='https://cloud.r-project.org')"

# Expose RStudio Server port
EXPOSE 8787

# Default working directory inside container
WORKDIR /home/rstudio

