#!/bin/bash
echo "SETTING LD LIBRARY PATH:"
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$(pwd)/dep/htslib-1.8/lib/:$(pwd)/dep/hdf5-1.10.2/build/hdf5/lib/
set LD_LIBRARY_PATH $LD_LIBRARY_PATH:$(pwd)/dep/htslib-1.8/lib/:$(pwd)/dep/hdf5-1.10.2/build/hdf5/lib/
echo $LD_LIBRARY_PATH