# Native Linux (arm64 or amd64) test image for tviblindi, with RStudio Server
# and vaevictis. See README.md in this directory for build/run instructions
# and why this exists.
FROM rocker/rstudio:4.5

RUN apt-get update && apt-get install -y --no-install-recommends \
      build-essential gfortran cmake \
      libcgal-dev libboost-dev libeigen3-dev libgmp-dev libmpfr-dev libglpk-dev \
      libcurl4-openssl-dev libssl-dev libxml2-dev libwebpmux3 \
      python3 python3-pip python3-venv python3-dev \
      wget ca-certificates \
    && rm -rf /var/lib/apt/lists/*

# Python venv for vaevictis, pinned to the commit documented in vaevictis/README.md
# ("pin to a specific commit rather than tracking master" - it's an explicitly unstable playground)
RUN python3 -m venv /opt/venv-tviblindi \
    && /opt/venv-tviblindi/bin/pip install --no-cache-dir --upgrade pip \
    && /opt/venv-tviblindi/bin/pip install --no-cache-dir tensorflow numpy scipy pandas scikit-learn \
    && /opt/venv-tviblindi/bin/pip install --no-cache-dir \
       "git+https://github.com/stuchly/vaevictis.git@50152ac30587a3b8f32ad67397526bcddf6a2f8e"

ENV RETICULATE_PYTHON=/opt/venv-tviblindi/bin/python

# tviblindi has no CGAL upper-version pin anymore (verified against 5.6-6.2;
# see tviblindi/README.md). Ubuntu Noble's apt libcgal-dev is 5.6-1build3,
# which needs BH pinned below 1.90.0 - newer BH bundles a Boost recent
# enough to have dropped boost::mpl::if_c (superseded by boost::mp11), which
# CGAL 5.x's headers (number_utils.h, Sqrt_extension_type.h, etc.) still
# reference. This pin is a CGAL-5.x-only concern - verified CGAL 6.x compiles
# fine even against BH 1.90+, since it no longer uses that Boost pattern.
RUN Rscript -e "install.packages('https://cran.r-project.org/src/contrib/Archive/BH/BH_1.87.0-1.tar.gz', repos=NULL, type='source')"

RUN Rscript -e "install.packages(c('Rcpp','RcppEigen','RcppArmadillo','shiny','shinyBS','shinybusy', \
      'shinycssloaders','shinyFiles','shinyWidgets','reticulate','tidyverse','gplots','ggdendro', \
      'dendextend','scattermore','cowplot','cluster','igraph','uwot','Matrix','jpeg','codetools', \
      'BiocManager'), repos='https://packagemanager.posit.co/cran/__linux__/noble/latest')" \
    && Rscript -e "BiocManager::install('flowCore', update=FALSE, ask=FALSE)"

# apt's libcgal-dev ships no compiled libCGAL.so (verified via `dpkg -L
# libcgal-dev` - header-only by design on this Ubuntu package). CGAL::Random
# and a few other normally-non-inline CGAL bits need the CGAL_HEADER_ONLY
# preprocessor macro defined to become inline instead of requiring that
# missing library. R's build system doesn't propagate a plain CPPFLAGS env
# var into the actual compile line (verified empirically), so inject it
# directly into the installed header instead.
RUN sed -i '15i #define CGAL_HEADER_ONLY 1' /usr/include/CGAL/config.h

COPY . /pkgsrc
RUN R CMD INSTALL /pkgsrc

COPY docker/linux-test/start_app.R /start_app.R

EXPOSE 8787 3838
