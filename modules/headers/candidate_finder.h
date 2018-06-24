//
// Created by Kishwar Shafin on 6/14/18.
//

#ifndef FRIDAY_CPP_CANDIDATE_FINDER_H
#define FRIDAY_CPP_CANDIDATE_FINDER_H

#include <iostream>
#include <sstream>
#include <set>
#include <stdio.h>
#include <string>
#include <vector>
#include <algorithm>
#include "vcf_handler.h"
#include "bam_handler.h"
#include "fasta_handler.h"

using namespace std;

#define GENOTYPE_HOM 0
#define GENOTYPE_HET 1
#define GENOTYPE_HOMALT 2


struct type_candidate_allele {
    long long pos;
    string allele;
    int candidate_type;
    int freq;
    int labeled_genotype;
    bool is_labeled;
    bool operator<(const type_candidate_allele& r) const {
        if(pos != r.pos) {
            return pos < r.pos;
        } else if(freq != r.freq){
            return freq > r.freq;
        } else if(candidate_type != r.candidate_type){
            return candidate_type < r.candidate_type;
        }
        return allele < r.allele;
    }
};

class candidate_finder {
    public:
        candidate_finder();
        map<long long, vector<type_candidate_allele> > find_candidates(string bam_file_path,
                                                                       string ref_file_path,
                                                                       string chromosome,
                                                                       long long start,
                                                                       long long stop,
                                                                       map<long long, int> &insert_length_map,
                                                                       string vcf_file_path=string(),
                                                                       bool label_candidates=false);
        int get_genotype(vector<int> genotype);
        ~candidate_finder();
};

#endif //FRIDAY_CPP_CANDIDATE_FINDER_H
