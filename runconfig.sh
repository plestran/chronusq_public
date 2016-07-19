#!/bin/sh

cmake \
 -DCMAKE_CXX_COMPILER=g++\
 -DCMAKE_C_COMPILER=gcc\
 -DCMAKE_PREFIX_PATH='/home/jjradler/local;/usr'\
 -DCMAKE_CXX_FLAGS='-w -O3'\
 -DBUILD_LIBINT=OFF\
 ..
