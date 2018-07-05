//
// Created by Kishwar Shafin on 7/5/18.
//

#ifndef FRIDAY_CPP_TRAIN_DATA_GENERATOR_H
#define FRIDAY_CPP_TRAIN_DATA_GENERATOR_H

#endif //FRIDAY_CPP_TRAIN_DATA_GENERATOR_H

#include <stdio.h>
#include <iostream>
#include <string.h>
#include <vector>
#include <algorithm>
#include "candidate_finder.h"
#include "bam_handler.h"
#include "bed_handler.h"
#include "image_generator.h"
#include "hdf5_handler.h"
#include "tqdm.h"
using namespace std;

#define ANSI_COLOR_RED     "\x1b[31m"
#define ANSI_COLOR_GREEN   "\x1b[32m"
#define ANSI_COLOR_YELLOW  "\x1b[33m"
#define ANSI_COLOR_BLUE    "\x1b[34m"
#define ANSI_COLOR_MAGENTA "\x1b[35m"
#define ANSI_COLOR_CYAN    "\x1b[36m"
#define ANSI_COLOR_RESET   "\x1b[0m"

class train_data_generator {
    public:
        train_data_generator(string bam_file_path,
                             string ref_file_path,
                             string vcf_file_path,
                             string confident_bed_path);
        void generate_labeled_images(string chromosome_name,
                                     long long start_pos,
                                     long long end_pos);
        void genome_level_processes();
        ~train_data_generator();
    private:
        string chromosome_name;
        string bam_file_path;
        string ref_file_path;
        string vcf_file_path;
        string confident_bed_path;
};
