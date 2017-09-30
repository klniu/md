FROM ubuntu:16.04
MAINTAINER Hugh Gao<hughgao01@gmail.com>

RUN mkdir -p ~/temp ~/tools/ \
    && cd ~/temp
    && apt-get update \
    && apt-get upgrade \
    && apt-get install libeigen2-dev libeigen3-dev python3 python3-dev cmake gcc git g++ swig2.0 \
    && wget https://github.com/openbabel/openbabel/archive/openbabel-2-4-1.tar.gz \
    && tar zxf openbabel-openbabel-2-4-1.tar.gz \
    && mkdir build \
    && cd build \
    && cmake ../openbabel-openbabel-2-4-1 -DCMAKE_INSTALL_PREFIX=~/tools/openbabel -DPYTHON_BINDINGS=ON -DRUN_SWIG=ON -DPYTHON_EXECUTABLE=/usr/bin/python3 \
    && make -j4 \
    && make install \
    && echo 'PATH=$PATH:~/tools/openbabel/bin:~/tools/md' >> ~/.bashrc \
    && echo 'export PYTHONPATH=${PYTHONPATH}:~/tools/openbabel/lib/python3.4/site-packages:~/tools/md' >> ~/.bashrc \
    && cd ~/tools \
    && git clone https://github.com/klniu/md.git \
