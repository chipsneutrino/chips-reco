#! /bin/bash

CURRENTDIR=$(pwd)

cd $CHIPSRECO
make clean
rm ./cmake_install.cmake
rm ./CMakeCache.txt
rm -r ./CMakeFiles
rm ./Makefile

cd $CURRENTDIR