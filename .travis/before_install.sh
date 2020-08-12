#!/bin/bash

if [ $TRAVIS_OS_NAME = 'osx' ]; then
	brew update;
	brew install hdf5;
	brew install eigen;
	brew install doxygen;
	brew install cmake;
else
    sudo apt-get libhdf5-dev;
    sudo apt-get libeigen3-dev;
    sudo apt-get doxygen;
    sudo apt-get lcov;
fi