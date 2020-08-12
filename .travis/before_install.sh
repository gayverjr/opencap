#!/bin/bash

if [ $TRAVIS_OS_NAME = 'osx' ]; then
	echo $CC
else
	sudo apt-get --force-yes install libhdf5-dev
	sudo apt-get --force-yes install libeigen3-dev
	sudo apt-get --force-yes install doxygen
	sudo apt-get --force-yes install lcov
fi