#!/bin/bash

mkdir -p include lib
rm -rf include/* lib/*

cd htslib
make
cp htslib/* ../include
cp libhts.a ../lib
cd -

make
