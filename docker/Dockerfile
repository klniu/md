FROM ubuntu:16.04
MAINTAINER Hugh Gao<hughgao01@gmail.com>

COPY ./sources.list /etc/apt/

ENV python_version 3.5
RUN mkdir -p ~/temp /opt/md \
    && cd ~/temp \
    && apt-get update \
    && apt-get -y upgrade \
    && apt-get -y install libeigen2-dev libeigen3-dev python3 python3-dev libcairo2-dev cmake gcc git g++ swig2.0 wget \
    && wget https://github.com/openbabel/openbabel/archive/openbabel-2-4-1.tar.gz \
    && tar zxf openbabel-2-4-1.tar.gz \
    && mkdir build \
    && cd build \
    && cmake ../openbabel-openbabel-2-4-1 -DCMAKE_INSTALL_PREFIX=/opt/openbabel -DPYTHON_BINDINGS=ON -DRUN_SWIG=ON -DPYTHON_EXECUTABLE=/usr/bin/python3 \
    && make -j4 \
    && make install \
    && git clone https://github.com/klniu/md.git /opt/md \
    && rm -r ~/temp \
    && apt-get -y clean && apt-get -y autoremove
# RUN py_version=`python3 --version|cut -c8-10` \
ENV PATH="/opt/openbabel/bin:/opt/md:${PATH}"
# ENV PYTHONPATH="/opt/openbabel/lib/python${py_version}/site-packages:/opt/md:${PYTHONPATH}"
ENV PYTHONPATH="/opt/openbabel/lib/python${python_version}/site-packages:/opt/md:${PYTHONPATH}"
