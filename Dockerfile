FROM rocker/tidyverse:4.4.1

RUN apt-get update && apt-get install -y libmagick++-dev

RUN R -e "install.packages(c('magick', 'devtools'))"
RUN R -e "BiocManager::install(c('SpatialExperiment', 'spaSim'), ask = FALSE)"
RUN R -e "devtools::install_github('WEHI-SODA-Hub/imcvis')"

CMD /bin/bash
