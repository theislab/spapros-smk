FROM bioconductor/bioconductor:devel

RUN apt-get update && \
    apt-get install -y libcurl4-openssl-dev libssl-dev libssh2-1-dev libxml2-dev zlib1g-dev

RUN R -e "install.packages('devtools')"
RUN R -e "devtools::install_github('MarioniLab/geneBasisR')"
RUN R -e "BiocManager::install('zellkonverter')"

# Install Python 3.7
RUN apt-get update
RUN apt-get install -y software-properties-common
RUN add-apt-repository ppa:deadsnakes/ppa
RUN apt-get install -y python3.7
RUN apt-get install -y python3-pip
RUN apt-get install -y python3.7-distutils

# Make Python 3.7 the default
RUN update-alternatives --install /usr/bin/python3 python3 /usr/bin/python3.7 1
RUN update-alternatives --install /usr/bin/python python /usr/bin/python3.7 1

RUN pip3 install scanpy