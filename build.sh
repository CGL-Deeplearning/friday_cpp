#!/bin/bash
make
echo "and the output is:"
export LD_LIBRARY_PATH=$(pwd)/dep/htslib-1.8/lib/:$(pwd)/dep/hdf5-1.10.2/build/hdf5/lib/:$LD_LIBRARY_PATH
set LD_LIBRARY_PATH $(pwd)/dep/htslib-1.8/lib/:$(pwd)/dep/hdf5-1.10.2/build/hdf5/lib/:$LD_LIBRARY_PATH
time ./friday.exe
echo $LD_LIBRARY_PATH