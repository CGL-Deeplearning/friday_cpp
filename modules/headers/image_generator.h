//
// Created by Kishwar Shafin on 6/23/18.
//

#ifndef FRIDAY_CPP_IMAGE_GENERATOR_H
#define FRIDAY_CPP_IMAGE_GENERATOR_H

#include <iostream>
#include <sstream>
#include <set>
#include <stdio.h>
#include <string>
#include <vector>
#include "candidate_finder.h"
#include "image_channels.h"
using namespace std;
#define CIGAR_OP_MATCH 0
#define CIGAR_OP_IN 1
#define CIGAR_OP_DEL 2

#define BASE_CHANNEL 0
#define BASE_QUAL_CHANNEL 1
#define MAP_QUAL_CHANNEL 2
#define STRAND_CHANNEL 3
#define MATCH_CHANNEL 4
#define CIGAR_CHANNEL 5
#define SUPPORT_CHANNEL 6

#define IMAGE_WIDTH 100
#define IMAGE_HEIGHT 50
#define TOTAL_CHANNELS 7


struct type_base_info {
    int base_quality;
    char base;
    bool is_match;
    bool operator<(const type_base_info& r) const {
        return base_quality < r.base_quality;
    }
};

class image_generator {
    public:
        image_generator(string chromosome_name,
                        string bam_file_path,
                        string ref_file_path,
                        type_candidate_allele candidate,
                        map<long long, int> insert_length_map);

        void set_left_right_genomic_position();
        void generate_candidate_image();
        void print_decoded_image();
        void set_reference_base(int row, int column, char base);
        void set_read_base(int row, int column,
                           char base, double base_qual,
                           double map_qual, bool is_rev,
                           bool is_match, bool is_support,
                           int cigar_op);
        ~image_generator();
        int image_array[IMAGE_HEIGHT][IMAGE_WIDTH][TOTAL_CHANNELS];
    private:
        string chromosome_name;
        string bam_file_path;
        string ref_file_path;
        type_candidate_allele candidate;
        map<long long, int> insert_length_map;
        int half_width;
        long long left_genomic_pos;
        long long right_genomic_pos;
        int total_left_bases;
        int total_right_bases;
};

#endif //FRIDAY_CPP_IMAGE_GENERATOR_H
