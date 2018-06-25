//
// Created by Kishwar Shafin on 6/23/18.
//

#include "../headers/image_generator.h"

void image_generator::set_left_right_genomic_position(){
    half_width = int(IMAGE_WIDTH / 2);

    // find the left and right positions where the image should start
    left_genomic_pos = candidate.pos - 1;
    total_left_bases = 0;
    while(total_left_bases < half_width){
        if(insert_length_map.find(left_genomic_pos) != insert_length_map.end()) {
            if (total_left_bases + insert_length_map[left_genomic_pos] <= half_width) {
                total_left_bases += (insert_length_map[left_genomic_pos]);
                left_genomic_pos -= 1;
            } else {
                left_genomic_pos += 1;
                break;
            }
        } else if(total_left_bases + 1 <= half_width) {
            total_left_bases += 1;
            left_genomic_pos -= 1;
        } else {
            left_genomic_pos += 1;
            break;
        }
    }

    right_genomic_pos = candidate.pos;
    total_right_bases = 0;
    while(total_right_bases < half_width){
        if(insert_length_map.find(right_genomic_pos) != insert_length_map.end()) {
            if (total_right_bases + insert_length_map[right_genomic_pos] <= half_width) {
                total_right_bases += (insert_length_map[right_genomic_pos]);
                right_genomic_pos += 1;
            } else {
                right_genomic_pos -= 1;
                break;
            }
        } else if(total_right_bases + 1 <= half_width) {
            total_right_bases += 1;
            right_genomic_pos += 1;
        } else {
            right_genomic_pos -= 1;
            break;
        }
    }
}

void image_generator::set_reference_base(int row, int column, char base) {
    image_channels get_channels;
    image_array[row][column][BASE_CHANNEL] = get_channels.get_base_color(base);
    image_array[row][column][BASE_QUAL_CHANNEL] = get_channels.get_base_quality_color(40);
    image_array[row][column][MAP_QUAL_CHANNEL] = get_channels.get_map_quality_color(60);
    image_array[row][column][STRAND_CHANNEL] = get_channels.get_strand_color(false);
    image_array[row][column][MATCH_CHANNEL] = get_channels.get_match_color(true);
    image_array[row][column][CIGAR_CHANNEL] = get_channels.get_cigar_color(CIGAR_OP_MATCH);
    image_array[row][column][SUPPORT_CHANNEL] = get_channels.get_support_color(true);
}

void image_generator::set_read_base(int row, int column,
                                    char base, double base_qual,
                                    double map_qual, bool is_rev,
                                    bool is_match, bool is_support,
                                    int cigar_op) {
    image_channels get_channels;
    image_array[row][column][BASE_CHANNEL] = get_channels.get_base_color(base);
    image_array[row][column][BASE_QUAL_CHANNEL] = get_channels.get_base_quality_color(base_qual);
    image_array[row][column][MAP_QUAL_CHANNEL] = get_channels.get_map_quality_color(map_qual);
    image_array[row][column][STRAND_CHANNEL] = get_channels.get_strand_color(is_rev);
    image_array[row][column][MATCH_CHANNEL] = get_channels.get_match_color(is_match);
    image_array[row][column][CIGAR_CHANNEL] = get_channels.get_cigar_color(cigar_op);
    image_array[row][column][SUPPORT_CHANNEL] = get_channels.get_support_color(is_support);
}

void image_generator::print_decoded_image() {
    //base channel decoding
    image_channels get_channels;
    for(int i=0; i<IMAGE_HEIGHT; i++) {
        for(int j=0; j<IMAGE_WIDTH; j++) {
            cout<<get_channels.get_decoded_base(image_array[i][j][BASE_CHANNEL]);
        }
        cout<<endl;
    }
}

void image_generator::generate_candidate_image() {
    set_left_right_genomic_position();
    BAM_handler bam_file(bam_file_path);

    // image array
    int current_column = half_width - total_left_bases;
    int current_row = 0;

    // First add the reference to the first row
    FASTA_handler fasta_file(ref_file_path);
    string reference_seq = fasta_file.get_reference_sequence(chromosome_name, left_genomic_pos, right_genomic_pos);

    //reference string
    for(int i=0; i<reference_seq.size(); i++) {
        this->set_reference_base(current_row, current_column, reference_seq[i]);
        current_column += 1;
        if(insert_length_map[this->left_genomic_pos+i] > 0){
            for(int l=0; l<insert_length_map[this->left_genomic_pos+i]; l++){
                this->set_reference_base(current_row, current_column, '*');
                current_column += 1;
            }
        }
    }
    current_row += 1;

    // get the id of the chromosome or sequence name
    const int tid = bam_name2id(bam_file.header, chromosome_name.c_str());

    // get the iterator
    hts_itr_t *iter  = sam_itr_queryi(bam_file.idx, tid, candidate.pos, candidate.pos+1);

    // initialize an alignment
    bam1_t* alignment = bam_init1();

    while(sam_itr_next(bam_file.hts_file, iter, alignment) >= 0) {
        //get read flags
        type_read_flags read_flags = bam_file.get_read_flags(alignment->core.flag);
        if(read_flags.is_supplementary || read_flags.is_qc_failed || read_flags.is_duplicate || read_flags.is_secondary
           || read_flags.is_unmapped){
            continue;
        }
        bool is_rev = read_flags.is_reverse;
        int map_qual = alignment->core.qual;
        // mapping quality
        if(map_qual <= 0){
            continue;
        }

        // get the base qualities and sequence bases
        uint8_t *seqi = bam_get_seq(alignment);
        uint8_t *qual = bam_get_qual(alignment);

        // get the cigar operations of the alignment
        uint32_t *cigar = bam_get_cigar(alignment);
        long long read_pos = alignment->core.pos;
        int read_index = 0;

        // initialize a map of read candidates
        map<long long, vector<type_candidate_allele> > read_candidate_map;
        map<long long, type_base_info> base_map;
        map<long long, string> insert_map;

        for(int k = 0; k < alignment->core.n_cigar; k++) {
            int cigar_op = bam_cigar_op(cigar[k]);
            int cigar_len = bam_cigar_oplen(cigar[k]);

            if(cigar_op == BAM_CMATCH){
                //match
                for(int i=0;i<cigar_len;i++) {
                    if(read_pos >= left_genomic_pos && read_pos <= right_genomic_pos) {
                        int ref_index = read_pos - left_genomic_pos;
                        bool is_match = true;
                        // get the read base
                        if(seq_nt16_str[bam_seqi(seqi, read_index)] != reference_seq[ref_index]) {
                            type_candidate_allele candidate;
                            candidate.pos = read_pos;
                            candidate.allele = seq_nt16_str[bam_seqi(seqi, read_index)];
                            candidate.candidate_type = SNP_TYPE;
                            read_candidate_map[read_pos].push_back(candidate);
                            is_match = false;
                        }
                        type_base_info base_info;
                        base_info.base = seq_nt16_str[bam_seqi(seqi, read_index)];
                        base_info.base_quality = int(qual[read_index]);
                        base_info.is_match = is_match;
                        base_map[read_pos] = base_info;
                    }
                    read_pos += 1;
                    read_index += 1;
                }
            } else if(cigar_op == BAM_CINS) {
                // insert
                long long anchor_position = read_pos - 1;
                if(anchor_position >= left_genomic_pos && anchor_position <= right_genomic_pos && read_index > 0) {
                    string insert_allele;
                    for (int i = -1; i < cigar_len; i++) {
                        insert_allele += seq_nt16_str[bam_seqi(seqi, read_index + i)];
                    }
                    type_candidate_allele candidate;
                    candidate.pos = anchor_position;
                    candidate.allele = insert_allele;
                    candidate.candidate_type = INSERT_TYPE;
                    read_candidate_map[anchor_position].push_back(candidate);

                    insert_map[anchor_position] = insert_allele;
                }
                read_index += cigar_len;
            } else if(cigar_op == BAM_CDEL) {
                // delete
                long long anchor_position = read_pos - 1;
                if(anchor_position >= left_genomic_pos && anchor_position + cigar_len <= right_genomic_pos && read_index > 0) {
                    // get the anchor allele from the read
                    string delete_allele;
                    delete_allele = seq_nt16_str[bam_seqi(seqi, read_index-1)];
                    int ref_index = read_pos - left_genomic_pos;
                    for (int i = 0; i < cigar_len; i++) {
                        delete_allele += reference_seq[ref_index + i];

                        type_base_info base_info;
                        base_info.base = '*';
                        base_info.base_quality = 20;
                        base_info.is_match = false;
                        base_map[read_pos + i] = base_info;
                    }
                    type_candidate_allele candidate;
                    candidate.pos = anchor_position;
                    candidate.allele = delete_allele;
                    candidate.candidate_type = DELETE_TYPE;
                    read_candidate_map[anchor_position].push_back(candidate);
                }
                read_pos += cigar_len;
            }else if(cigar_op == BAM_CREF_SKIP){
                read_pos += cigar_len;
            }else if(cigar_op == BAM_CSOFT_CLIP){
                read_index += cigar_len;
            }else if(cigar_op == BAM_CHARD_CLIP){
                // hard clip
            }else if(cigar_op == BAM_CPAD){
                // padding
            }else if(cigar_op == BAM_CEQUAL){
                // is a match
                read_pos += cigar_len;
                read_index += cigar_len;
            }else if(cigar_op == BAM_CDIFF){
                // is a mismatch
                read_pos += cigar_len;
                read_index += cigar_len;
            }else if(cigar_op == BAM_CBACK){
                cerr<<"BAM CONTAINS BAM_CBACK WHICH WE DON'T SUPPORT YET"<<endl;
            }
        }

        bool is_support = false;
        vector<type_candidate_allele> read_candidates = read_candidate_map[candidate.pos];
        // this will help to find the support channels
        for(int i=0; i<read_candidates.size(); i++) {
            if(candidate.allele ==  read_candidates[i].allele &&
               candidate.candidate_type == read_candidates[i].candidate_type) {
                is_support = true;
                break;
            }
        }

        // set everything to the image row
        uint32_t len = alignment->core.l_qseq;
        long long read_alignment_pos = alignment->core.pos;
        long long read_alignment_end_pos = bam_endpos(alignment);
        current_column = half_width - total_left_bases;
        if(read_alignment_pos > left_genomic_pos)
            current_column = current_column + (read_alignment_pos-left_genomic_pos);

        for(long long i=read_alignment_pos; i<read_alignment_end_pos; i++) {
            if(i < left_genomic_pos) {
                continue;
            }
            if( i >= right_genomic_pos) {
                break;
            }
            if(base_map.find(i) != base_map.end()) {
                char base = base_map[i].base;
                int base_qual = base_map[i].base_quality;
                bool is_match = base_map[i].is_match;
                int cigar_op = CIGAR_OP_MATCH;
                if(base == '*') {
                    cigar_op = CIGAR_OP_DEL;
                }
                if(current_column < IMAGE_WIDTH) {
                    set_read_base(current_row, current_column, base, base_qual, map_qual, is_rev, is_match, is_support,
                                  cigar_op);
                    current_column += 1;
                } else {
                    break;
                }
            }

            if(insert_length_map.find(i) != insert_length_map.end()) {
                int inserted_bases = 0;
                if(insert_map.find(i) != insert_map.end()) {
                    string insert_allele = insert_map[i];
                    for(int i=1; i<insert_allele.size(); i++) {
                        char base = insert_allele[i];
                        int base_qual = 20;
                        bool is_match = false;
                        int cigar_op = CIGAR_OP_IN;
                        if(current_column < IMAGE_WIDTH) {
                            set_read_base(current_row, current_column, base, base_qual, map_qual, is_rev, is_match,
                                          is_support, cigar_op);
                            current_column += 1;
                        } else {
                            break;
                        }
                    }
                    inserted_bases = insert_allele.size();
                }
                int extra_in_bases = insert_length_map[i] - inserted_bases;

                for(int i=0; i<extra_in_bases; i++) {
                    char base = '*';
                    int base_qual = 20;
                    bool is_match = false;
                    int cigar_op = CIGAR_OP_IN;
                    if(current_column < IMAGE_WIDTH) {
                        set_read_base(current_row, current_column, base, base_qual, map_qual, is_rev, is_match,
                                      is_support, cigar_op);
                        current_column += 1;
                    } else {
                        break;
                    }
                }
            }
        }
        current_row += 1;
        if(current_row >= IMAGE_HEIGHT)
            break;

    }

    hts_itr_destroy(iter);
    bam_destroy1(alignment);
    print_decoded_image();
}

image_generator::image_generator(string chromosome_name, string bam_file_path, string ref_file_path,
                                 type_candidate_allele candidate, map<long long, int> insert_length_map) {
    this->chromosome_name = chromosome_name;
    this->bam_file_path = bam_file_path;
    this->ref_file_path = ref_file_path;
    this->candidate = candidate;
    this->insert_length_map = insert_length_map;
    memset(this->image_array, 0, sizeof(this->image_array));
}

image_generator::~image_generator() {
}