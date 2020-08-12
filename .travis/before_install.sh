#!/bin/bash

if [ $TRAVIS_OS_NAME = 'osx' ]; then
	brew update
	brew install hdf5
	brew install eigen
	brew install doxygen
	brew install cmake
else
	apt-get libhdf5-dev
	apt-get libeigen3-dev
	apt-get doxygen
	apt-get lcov
fi