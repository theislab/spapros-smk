# Use the Rocker r-ver image with R version 4.0.3
FROM rocker/r-ver:4.0.3

# Update the package list
RUN apt-get update
RUN apt-get install -y libpng-dev

# Install required packages
RUN R -e "install.packages('reticulate', repos='https://cran.rstudio.com/')"
RUN R -e "install.packages('pracma', repos='https://cran.rstudio.com/')"

# Install Python 3.7
RUN apt-get install -y software-properties-common
RUN add-apt-repository ppa:deadsnakes/ppa
RUN apt-get install -y python3.7
RUN apt-get install -y python3-pip
RUN apt-get install -y python3.7-distutils
RUN apt-get install -y libpython3.7-dev

# Make Python 3.7 the default
RUN update-alternatives --install /usr/bin/python3 python3 /usr/bin/python3.7 1
RUN update-alternatives --install /usr/bin/python python /usr/bin/python3.7 1

RUN pip3 install scanpy

# Install R packages that I forgot first
RUN R -e "install.packages('anndata')"

