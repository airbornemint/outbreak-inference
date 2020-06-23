############################################################
# R with packages and libraries we need
FROM rocker/r-ver AS r
LABEL maintainer="Ben Artin <ben@artins.org>"

### Setup apt packages needed to build the image
ARG DEBIAN_FRONTEND=noninteractive
RUN apt-get update && apt-get install --yes --no-install-recommends moreutils > /dev/null 2>&1
# SHELL [ "/usr/bin/chronic", "/bin/bash", "-c" ]

# Most of these are needed to build a dependency of the package, which means they 
# are needed for to build dependencies of the paper too
RUN apt-get update && apt-get install --yes --no-install-recommends \
	libcurl4-gnutls-dev \
	gnutls-dev \
	libssh2-1-dev \
	libxml2-dev \
	zlib1g-dev \
	libpng-dev \
	libgit2-dev \
	libssl-dev

# Pandoc is needed for vignettes	
RUN apt-get update && apt-get install --yes --no-install-recommends \
	pandoc \
	pandoc-citeproc \
	libxt

############################################################
# R with dependencies for building the package
FROM r AS package-tools

WORKDIR /package

# Deps
COPY .Rprofile .Rprofile
COPY renv renv
# Devtools before deps to avoid rebuilding devtools if deps change
RUN Rscript -e "renv::install('devtools')"
COPY renv.lock renv.lock
RUN Rscript -e "renv::restore()"
COPY . .

SHELL ["/bin/bash", "-c"]
ENTRYPOINT ["Rscript", "-e"]

############################################################
# Tex environment we use to build the paper (it includes R because of knitr)
FROM r AS tex
LABEL maintainer="Ben Artin <ben@artins.org>"

RUN apt-get update && apt-get install --yes --no-install-recommends \
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
	texlive-pictures \
	latexmk \
	poppler-utils \
	imagemagick

RUN tlmgr init-usertree
RUN tlmgr --usermode option repository http://ftp.math.utah.edu/pub/tex/historic/systems/texlive/2017/tlnet-final/

RUN luaotfload-tool -u

WORKDIR /tex
ENTRYPOINT ["/bin/bash", "-c"]

############################################################
# Commands for knitting the paper source
FROM tex AS paper-tools

# Install dependencies
WORKDIR /tex

COPY Paper/.Rprofile .Rprofile
COPY Paper/renv renv
COPY Paper/renv.lock renv.lock
RUN Rscript -e "renv::restore()"

COPY Paper .
COPY vignettes/seasonal.csv /vignettes/seasonal.csv
