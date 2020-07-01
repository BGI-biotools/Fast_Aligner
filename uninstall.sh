#!/bin/bash

rm -rf include lib

cd htslib
make clean
cd -

make clean
