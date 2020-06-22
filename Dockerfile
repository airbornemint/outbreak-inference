############################################################
# R with packages and libraries we need
FROM rocker/r-ver AS r
LABEL maintainer="Ben Artin <ben@artins.org>"

### Setup apt packages needed to build the image
ARG DEBIAN_FRONTEND=noninteractive
RUN apt-get update && apt-get install --yes --no-install-recommends moreutils > /dev/null 2>&1
SHELL [ "/usr/bin/chronic", "/bin/bash", "-c" ]

### Setup R packages
RUN apt-get update && apt-get install --yes --no-install-recommends \
	libcurl4-gnutls-dev \
	gnutls-dev \
	libssh2-1-dev \
	libxml2-dev \
	zlib1g-dev \
	libpng-dev \
	libgit2-dev \
	libssl-dev

# install2.r needs these
RUN Rscript -e "install.packages(c('docopt', 'remotes'))"

# These come from CRAN
RUN install2.r \
	zipcode \
	dplyr \
	knitr \
	reshape \
	mgcv \
	data.table \
	tikzDevice \
	sp \
	mapproj \
	ggplot2 \
	ggstance \
	gridExtra \
	devtools \
	import \
	doParallel \
	kableExtra \
	rmarkdown \
	plotrix

WORKDIR /package
SHELL ["/bin/bash", "-c"]
ENTRYPOINT ["/bin/bash", "-c"]
CMD ["/bin/bash"]

############################################################
# Commands for building the R package
FROM r AS build-package

COPY . /package

SHELL [ "/usr/bin/chronic", "/bin/bash", "-c" ]
ENTRYPOINT Rscript -e "devtools::install_dev_deps(upgrade=FALSE); devtools::build(pkg='/package', path=Sys.getenv('R_BUILD_DIR'), vignettes=as.logical(Sys.getenv('R_BUILD_VIGNETTES')), binary=FALSE)"

############################################################
# Commands for checking the R package
FROM r AS check-package

# For dev deps
COPY DESCRIPTION . 

# SHELL [ "/usr/bin/chronic", "/bin/bash", "-c" ]
ENTRYPOINT Rscript -e "devtools::install_dev_deps(upgrade=FALSE); devtools::check_built(path=Sys.getenv('R_PACKAGE_ARCHIVE'))"

############################################################
# Tex environment we use (it includes R because knitr needs both)
FROM r AS tex
LABEL maintainer="Ben Artin <ben@artins.org>"

SHELL [ "/usr/bin/chronic", "/bin/bash", "-c" ]
RUN apt-get update && apt-get install --yes --no-install-recommends \
	pandoc \
	pandoc-citeproc \
	qpdf \
	wget \
	xzdec \
	lmodern \
	texlive \
	texlive-binaries \
	texlive-luatex \
	texlive-lang-cyrillic \
	texlive-latex-extra \
	texlive-bibtex-extra \
	texlive-fonts-extra \
	texlive-pictures

RUN tlmgr init-usertree
RUN tlmgr --usermode option repository http://ftp.math.utah.edu/pub/tex/historic/systems/texlive/2017/tlnet-final/

WORKDIR /tex
SHELL [ "/bin/bash", "-c" ]
ENTRYPOINT /bin/bash

############################################################
# Commands for knitting the paper source
# Knitr requires latex for masuring figure elements
FROM tex AS knitr-paper

SHELL [ "/usr/bin/chronic", "/bin/bash", "-c" ]

# Install dependencies
WORKDIR /tex

COPY Paper/.Rprofile .Rprofile
COPY Paper/renv renv
COPY Paper/renv.lock renv.lock
RUN Rscript -e "renv::restore()"

COPY Paper /tex
COPY vignettes/seasonal.csv /vignettes/seasonal.csv

SHELL [ "/bin/bash", "-c" ]
ENTRYPOINT Rscript -e "install.packages(Sys.getenv('R_PACKAGE_ARCHIVE')); options(pspline.paper.validation.run=as.logical(Sys.getenv('KNITR_RUN_VALIDATION'))); options(pspline.paper.output=Sys.getenv('KNITR_OUTPUT_DIR')); knitr::knit(Sys.getenv('KNITR_INPUT_FILE'), output=Sys.getenv('KNITR_OUTPUT_FILE'))"

############################################################
# Commands for converting figures to PDF
FROM knitr-paper AS pdflatex-figures

SHELL [ "/usr/bin/chronic", "/bin/bash", "-c" ]
ENTRYPOINT cd "${KNITR_OUTPUT_DIR}/figures" ; find . -name "*.tex" -exec pdflatex {} \;

############################################################
# Commands for generating paper PDF
FROM knitr-paper AS pdflatex-paper

RUN apt-get update && apt-get install --yes --no-install-recommends \
	latexmk

ENTRYPOINT latexmk -pdflua "${LATEX_FILE}"