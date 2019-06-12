FROM rocker/r-ver
LABEL maintainer="Ben Artin <ben@artins.org>"

### Setup apt packages needed to build the image
ARG DEBIAN_FRONTEND=noninteractive
RUN apt-get update -y -qq

# These are only needed to build the image (and will not be present in the final image)
# ARG BUILD_APT_DEPENDENCIES=""
# RUN apt-get install -y -qq --no-install-recommends ${BUILD_APT_DEPENDENCIES}
# RUN apt-mark auto ${BUILD_APT_DEPENDENCIES}

### Setup R packages
ARG R_APT_DEPENDENCIES="libcurl4-gnutls-dev gnutls-dev libssh2-1-dev libxml2-dev zlib1g-dev libpng-dev libgit2-dev libssl-dev pandoc pandoc-citeproc qpdf"
RUN apt-get install -y -qq --no-install-recommends ${R_APT_DEPENDENCIES}

### Setup misc other packages needed to run our build
ARG RUN_APT_DEPENDENCIES="tzdata locales gpg"
RUN apt-get install -y -qq --no-install-recommends ${RUN_APT_DEPENDENCIES}

# Finish up the configuration
RUN echo tzdata tzdata/Areas select America >> /tmp/debconf.txt
RUN echo tzdata tzdata/Zones/America select New_York >> /tmp/debconf.txt
RUN debconf-set-selections /tmp/debconf.txt
RUN rm /etc/timezone
RUN rm /etc/localtime
RUN dpkg-reconfigure --frontend noninteractive tzdata
RUN update-locale --reset LANG=en_US.utf8 LANGUAGE=en_US LC_ALL=en_US.utf8 LC_CTYPE=en_US.utf8

# Config R before we go installing packages
RUN echo "options(Ncpus = $(nproc --all))" >> "/usr/local/lib/R/etc/Rprofile.site"

# install2.r needs these
RUN Rscript -e "install.packages(c('docopt', 'remotes'))"

# These come from CRAN
ARG BUILD_CRAN_DEPENDENCIES="zipcode dplyr knitr reshape mgcv data.table tikzDevice sp mapproj ggplot2 ggstance gridExtra devtools import doParallel kableExtra rmarkdown"
RUN /usr/local/lib/R/site-library/littler/examples/install2.r ${BUILD_CRAN_DEPENDENCIES}

# These come from GitHub
# ARG BUILD_R_GITHUB_DEPENDENCIES="airbornemint/outbreak-inference"
# RUN /usr/local/lib/R/site-library/littler/examples/installGithub.r ${BUILD_R_GITHUB_DEPENDENCIES}

### Setup TeX packages
ARG TEX_APT_DEPENDENCIES="wget xzdec lmodern texlive texlive-binaries texlive-luatex texlive-lang-cyrillic texlive-latex-extra texlive-bibtex-extra texlive-fonts-extra texlive-pictures"
RUN apt-get install -y -qq --no-install-recommends ${TEX_APT_DEPENDENCIES}

RUN tlmgr init-usertree
RUN tlmgr --usermode option repository http://ftp.math.utah.edu/pub/tex/historic/systems/texlive/2017/tlnet-final/

# ARG BUILD_TEX_DEPENDENCIES="xcharter"
# RUN tlmgr --usermode install ${BUILD_TEX_DEPENDENCIES}

### Cleanup to reduce size of the image
RUN apt-get autoremove -y
RUN apt-get clean
RUN rm -rf /tmp/*

ENTRYPOINT /bin/bash