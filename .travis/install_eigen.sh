#!/bin/bash

git clone https://gitlab.com/libeigen/eigen.git
git checkout tags/3.3.7
cd eigen
cmake .
make install
cd ..