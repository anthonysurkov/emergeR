# Dockerfile
FROM rocker/r2u:latest

# System + R deps via apt
RUN apt-get update && apt-get install -y --no-install-recommends \
	r-cran-dplyr r-cran-purrr r-cran-tidyr r-cran-tibble r-cran-magrittr \
    r-cran-tidyselect r-cran-binom r-cran-ade4 r-cran-rtsne r-cran-dbscan \
    r-cran-ggplot2 r-cran-ggprism r-cran-ggseqlogo r-cran-scales r-cran-rlang \
	r-bioc-biostrings r-bioc-iranges r-bioc-s4vectors \
    && apt-get clean && rm -rf /var/lib/apt/lists/*

RUN R -e "install.packages('devtools', repos = 'https://cloud.r-project.org')"

# Copy pkg src
WORKDIR /pkg
COPY . /pkg

CMD ["R"]
