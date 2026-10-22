# Linux R check image for scripts/check-container-local.sh.
# Build context is this directory plus a copy of DESCRIPTION, so the dependency
# layer is cached until the dependency list changes.
ARG R_TAG=4.5
FROM rocker/r-ver:${R_TAG}

# Vignettes use the quarto engine; install_quarto.sh also provides pandoc.
RUN /rocker_scripts/install_quarto.sh && /rocker_scripts/install_pandoc.sh

RUN apt-get update \
 && apt-get install -y --no-install-recommends \
      libcurl4-openssl-dev libssl-dev libxml2-dev libglpk-dev libnlopt-dev \
      libfontconfig1-dev libfreetype6-dev libharfbuzz-dev libfribidi-dev \
      libpng-dev libtiff5-dev libjpeg-dev qpdf \
 && rm -rf /var/lib/apt/lists/*

COPY DESCRIPTION /tmp/DESCRIPTION
RUN Rscript -e 'install.packages(c("remotes", "rcmdcheck"))' \
 && Rscript -e 'remotes::install_deps("/tmp", dependencies = TRUE)'

WORKDIR /pkg
