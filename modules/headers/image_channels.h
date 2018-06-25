//
// Created by Kishwar Shafin on 6/23/18.
//

#ifndef FRIDAY_CPP_IMAGE_CHANNELS_H
#define FRIDAY_CPP_IMAGE_CHANNELS_H

#include <algorithm>
using namespace std;
#define MAX_COLOR_VALUE 254.0
#define BASE_QUALITY_CAP 40.0
#define MAP_QUALITY_CAP 60.0
#define MAP_QUALITY_FILTER 5.0
#define MIN_DELETE_QUALITY 20.0

class image_channels{
    public:
    char get_decoded_base(double base_value);
    char get_decoded_support(double support_value);

    double get_base_color(char base);
    double get_base_quality_color(int qual);
    double get_map_quality_color(int qual);
    double get_strand_color(bool is_rev);
    double get_match_color(bool is_match);
    double get_cigar_color(int cigar_op);
    double get_support_color(bool is_supported);
};
#endif //FRIDAY_CPP_IMAGE_CHANNELS_H
