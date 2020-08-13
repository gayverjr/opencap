#!/bin/bash
if [ $TRAVIS_OS_NAME = 'linux' ]; then
	sudo apt-get --force-yes install libhdf5-dev
	sudo apt-get --force-yes install libeigen3-dev
	sudo apt-get --force-yes install doxygen
	sudo apt-get --force-yes install lcov
fi