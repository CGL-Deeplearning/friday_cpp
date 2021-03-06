//
// Created by Kishwar Shafin on 6/12/18.
//

#ifndef FRIDAY_CPP_BAM_HANDLER_H
#define FRIDAY_CPP_BAM_HANDLER_H

#include <iostream>
#include <sstream>
#include <set>
#include <stdio.h>
#include <string>
#include <vector>
#include <map>
#include <htslib/sam.h>

using namespace std;

#define SNP_TYPE 1
#define INSERT_TYPE 2
#define DELETE_TYPE 3

typedef struct{
    int cigar_op;
    int cigar_len;
} type_cigar;

typedef struct type_read_flags{
    bool is_paired;
    bool is_proper_pair;
    bool is_unmapped;
    bool is_mate_unmapped;
    bool is_reverse;
    bool is_mate_is_reverse;
    bool is_read1;
    bool is_read2;
    bool is_secondary;
    bool is_qc_failed;
    bool is_duplicate;
    bool is_supplementary;
    type_read_flags(){
        is_paired= 0;
        is_proper_pair= 0;
        is_unmapped= 0;
        is_mate_unmapped= 0;
        is_reverse= 0;
        is_mate_is_reverse= 0;
        is_read1= 0;
        is_read2= 0;
        is_secondary= 0;
        is_qc_failed= 0;
        is_duplicate= 0;
        is_supplementary= 0;
    }
} type_read_flags;

typedef struct{
    long long pos;
    string read_id;
    type_read_flags flags;
	string sequence;
	vector <type_cigar> cigar_sequence;
	int mapping_quality;
	vector <int> base_qualities;
} type_read;

typedef struct{
    string sequence_name;
    int sequence_length;
} type_sequence;


class BAM_handler {
    public:
        htsFile* hts_file;
        bam_hdr_t* header;
        hts_idx_t* idx;

        BAM_handler(string path);

        vector<type_read> get_reads(string region, long long start, long long stop);
        vector<string> get_chromosome_sequence_names();
        vector<type_sequence> get_chromosome_sequence_names_with_length();
        set<string> get_sample_names();
        type_read_flags get_read_flags(int flag);

    	~BAM_handler();

};

#endif // BAM_HANDLER_H
