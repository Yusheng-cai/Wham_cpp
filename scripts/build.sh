#!/bin/bash

export CC=gcc-11
export CXX=g++-11

# specify the build directory
build_type=RELEASE
build_dir=$PWD/${build_type}/
install_dir=$HOME/programs/Wham
fftw_dir=${HOME}/programs/fftw-3.3.10

# remove build_dir if it already exists
if [ -d $build_dir ] 
    then 
    rm -r $build_dir
fi

if [ -d $install_dir ]
    then 
    rm -rf $install_dir
fi

# make the build directory
mkdir -p $build_dir

# configure the build with cmake
cd $build_dir
cmake .. \
	-DCMAKE_BUILD_TYPE=${build_type} \
	-DCMAKE_INSTALL_PREFIX=${install_dir} \
	-DFFTW3_DIR=${fftw_dir}

# make with 8 threads
make -j 8 

make test
make install
