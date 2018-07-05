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
using namespace std;

class hdf5_handler{
    public:
        hdf5_handler(const string file_name);
        void save_img_to_hdf5(const int candidate_image[IMAGE_HEIGHT][IMAGE_WIDTH][TOTAL_CHANNELS]);
    void save_allele_info(string chromosome_name, string position,
                                        string allele1, string allele1_type,
                                        string allele2, string allele2_type);
        ~hdf5_handler();
    private:
        string file_name;

};

#endif //FRIDAY_CPP_HDF5_HANDLER_H
