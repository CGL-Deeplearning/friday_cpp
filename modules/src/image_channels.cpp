//
// Created by Kishwar Shafin on 6/24/18.
//

#include "../headers/image_channels.h"

char image_channels::get_decoded_base(double base_value){
    if(base_value >= 245.0 && base_value <= 260.0)
        return 'A';
    if(base_value >= 95.0 && base_value <= 105.0)
        return 'C';
    if(base_value >= 175.0 && base_value <= 185.0)
        return 'G';
    if(base_value >= 25.0 && base_value <= 35.0)
        return 'T';
    if(base_value >= 4.0 && base_value <= 10.0)
        return '*';
    return ' ';
}

double image_channels::get_base_color(char base){
    if(base == 'A' || base == 'a')
        return 254.0;
    if(base == 'C' || base == 'c')
        return 100.0;
    if(base == 'G' || base=='g')
        return 180.0;
    if(base == 'T' || base == 't')
        return 30.0;
    if(base == '*' || base == 'N' )
        return 5.0;
    return 1.0;
}

double image_channels::get_base_quality_color(int base_quality){
    if(base_quality > BASE_QUALITY_CAP) base_quality = BASE_QUALITY_CAP;

    return (MAX_COLOR_VALUE * base_quality) / BASE_QUALITY_CAP;
}

double image_channels::get_map_quality_color(int map_quality){
    if(map_quality > MAP_QUALITY_CAP) map_quality = MAP_QUALITY_CAP;
    return (MAX_COLOR_VALUE * map_quality) / MAP_QUALITY_CAP;
}

double image_channels::get_strand_color(bool is_rev){
    if(is_rev)
        return MAX_COLOR_VALUE;
    return MAX_COLOR_VALUE * 0.5;
}

double image_channels::get_match_color(bool is_match){
    if (is_match)
        return MAX_COLOR_VALUE * 0.2;
    else
        return MAX_COLOR_VALUE * 1.0;
}

double image_channels::get_cigar_color(int cigar_op){
    if(cigar_op == 2)
        return MAX_COLOR_VALUE;
    if(cigar_op == 1)
        return MAX_COLOR_VALUE * 0.6;

    return MAX_COLOR_VALUE * 0.3;
}

double image_channels::get_support_color(bool is_supported){
    if (is_supported)
        return MAX_COLOR_VALUE * 0.2;
    else
        return MAX_COLOR_VALUE * 1.0;
}