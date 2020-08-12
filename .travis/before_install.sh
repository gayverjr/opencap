#!/bin/bash

if [ $TRAVIS_OS_NAME = 'osx' ]; then
	brew update
	brew install hdf5
	brew install eigen
	brew install gcc@9
else
    if [ "$TRAVIS_PYTHON_VERSION" = "3.6" ]
    then
    	sudo apt-get --force-yes install python3.8-venv
    	sudo apt-get --force-yes install python3.8-dev
    	sudo apt-get --force-yes install python3.8-dev
    elif [ "$TRAVIS_PYTHON_VERSION" = "3.7" ]
    then
    	sudo apt-get --force-yes install python3.7-venv
    	sudo apt-get --force-yes install python3.7-dev
    	sudo apt-get --force-yes install python3.7-dev
    elif [ "$TRAVIS_PYTHON_VERSION" = "3.8" ]
    then
    	sudo apt-get --force-yes install python3.8-venv
    	sudo apt-get --force-yes install python3.8-dev
    	sudo apt-get --force-yes install python3.8-dev
    fi
	sudo apt-get --force-yes install libhdf5-dev
	sudo apt-get --force-yes install libeigen3-dev
	sudo apt-get --force-yes install doxygen
	sudo apt-get --force-yes install lcov

fi