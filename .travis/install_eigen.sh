#!/bin/bash

git clone https://gitlab.com/libeigen/eigen.git
cd eigen
git checkout tags/3.3.7
mkdir build
cd build
cmake ..
make install
cd ../../