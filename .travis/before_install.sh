#!/bin/bash

if [ $TRAVIS_OS_NAME = 'osx' ]; then
	brew update
	brew install hdf5
	brew install eigen
	brew install doxygen
	brew install gcc@9
else
	sudo apt-get --force-yes install libhdf5-dev
	sudo apt-get --force-yes install libeigen3-dev
	sudo apt-get --force-yes install doxygen
	sudo apt-get --force-yes install lcov
	sudo apt-get --force-yes install python3.8-venv
    sudo apt-get --force-yes install python3.8-dev
    sudo apt-get --force-yes install python3.8-dev
fi