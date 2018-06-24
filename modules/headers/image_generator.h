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
#include "H5Cpp.h"
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
#define TOTAL_CHANNELS 7
class image_generator {
    public:
        image_generator(string bam_file_path,
                        string ref_file_path,
                        string vcf_file_path=string(),
                        bool train_mode=false);

        void parse_candidates(string chromosome_name, long long start, long long stop);
        void generate_candidate_image(string chromosome_name,
                                      type_candidate_allele candidate,
                                      map<long long, int> &insert_length_map,
                                      const int image_width=300,
                                      const int image_height=100);
        ~image_generator();
    private:
        string bam_file_path;
        string ref_file_path;
        string vcf_file_path;
        bool is_train_mode;
};

#endif //FRIDAY_CPP_IMAGE_GENERATOR_H
