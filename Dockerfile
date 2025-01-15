FROM rocker/tidyverse:4.4.1

ENV XDG_CACHE_HOME=/tmp/.cache

RUN apt-get update && \
    apt-get install -y libmagick++-dev && \
    apt-get install -y curl && \
    curl -LO https://quarto.org/download/latest/quarto-linux-amd64.deb && \
    dpkg -i quarto-linux-amd64.deb && \
    rm quarto-linux-amd64.deb

RUN R -e "install.packages(c('magick', 'devtools', 'quarto'))"
RUN R -e "BiocManager::install(c('SpatialExperiment', 'spaSim'), ask = FALSE)"
RUN R -e "devtools::install_github('WEHI-SODA-Hub/spatialVis')"

CMD /bin/bash
