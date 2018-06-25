//
// Created by Kishwar Shafin on 6/25/18.
//

#ifndef FRIDAY_CPP_HDF5_HANDLER_H
#define FRIDAY_CPP_HDF5_HANDLER_H

#include <iostream>
#include <sstream>
#include <set>
#include <stdio.h>
#include <string>
#include <vector>
#include "candidate_finder.h"
#include "image_generator.h"
#include "hdf5.h"
#include "H5Cpp.h"
using namespace H5;
using namespace std;

class hdf5_handler{
    public:
        void save_img_to_hdf5(image_generator candidate_image);
};

#endif //FRIDAY_CPP_HDF5_HANDLER_H
