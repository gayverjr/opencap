#!/bin/bash
wget https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.10/hdf5-1.10.5/src/hdf5-1.10.5.tar.gz
tar -zxvf hdf5-1.10.5.tar.gz
cd hdf5-1.10.5
./configure
make install