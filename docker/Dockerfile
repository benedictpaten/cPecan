FROM ubuntu:14.04
MAINTAINER Trevor Pesout, tpesout@ucsc.edu

# update and install dependencies (not sure if all these are needed, but it works!)
RUN apt-get update && apt-get -y install git make gcc bzip2 lzma-dev zlib1g-dev libcurl4-openssl-dev libcrypto++-dev \
    libpthread-stubs0-dev g++ wget libbz2-dev liblzma-dev autoconf unzip php5-mcrypt python-pip python-dev
RUN pip install numpy pysam

# get sonLib
WORKDIR /opt
ENV SONLIB_COMMIT 265d6c510cff90412d88e5ab9945818cfdff6e41
RUN wget https://github.com/benedictpaten/sonLib/archive/$SONLIB_COMMIT.zip && unzip $SONLIB_COMMIT.zip \
    && ls && cd /opt/sonLib-$SONLIB_COMMIT && make && cd /opt && mv /opt/sonLib-$SONLIB_COMMIT /opt/sonLib \
    && rm -rf /opt/*$SONLIB_COMMIT*

# get cPecan
WORKDIR /opt
ADD tmp.cPecan /opt/cPecan

# build cPecan
WORKDIR /opt/cPecan
RUN make clean ; make

# setup entrypoint
COPY wrapper.sh /opt/
ENTRYPOINT ["sh", "/opt/wrapper.sh"]
