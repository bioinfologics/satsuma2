#! /bin/bash
# Bash Script that builds project

mkdir build
cd build
mkdir product
cmake -DCMAKE_INSTALL_PREFIX=product .. ${CMAKE_OPTIONS}
make all -j8
make install

tar cz product > satsuma2-${TRAVIS_OS_NAME}.tar.gz


