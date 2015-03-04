#!/bin/bash

curl -OL http://ftpmirror.gnu.org/autoconf/autoconf-2.69.tar.gz
curl -OL http://ftpmirror.gnu.org/automake/automake-1.14.tar.gz
curl -OL http://ftpmirror.gnu.org/libtool/libtool-2.4.2.tar.gz

tar xvf autoconf-2.69.tar.gz
tar xvf automake-1.14.tar.gz
tar xvf libtool-2.4.2.tar.gz

cd autoconf-2.69
./configure && make && sudo make install
cd ..
cd automake-1.14
./configure FC=gfortran F77=gfortran && make && sudo make install
cd ..
cd libtool-2.4.2
./configure FC=gfortran F77=gfortran && make && sudo make install
cd ..
