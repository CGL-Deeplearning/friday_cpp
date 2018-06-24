//
// Created by Kishwar Shafin on 6/23/18.
//

#include "../headers/image_generator.h"

void image_generator::generate_candidate_image(string chromosome_name,
                                               type_candidate_allele candidate,
                                               map<long long, int> &insert_length_map,
                                               const int image_width,
                                               const int image_height) {
    BAM_handler bam_file(this->bam_file_path);
    int half_width = int(image_width / 2);

    // find the left and right positions where the image should start
    int left_pos = candidate.pos - 1;
    int left_pixels = 0;
    while(left_pixels < half_width){
        if(left_pixels + insert_length_map[left_pos] + 1 <half_width) {
            left_pixels += ( insert_length_map[left_pos] + 1);
            left_pos -= 1;
        }
        else break;
    }
    int right_pos = candidate.pos;
    int right_pixels = 0;
    while(right_pixels < half_width){
        if(right_pixels + insert_length_map[right_pos] + 1 <= half_width) {
            right_pixels += ( insert_length_map[right_pos] + 1);
            right_pos += 1;
        }
        else break;
    }
    // image array
    uint8_t image_array[image_height][image_width][7];
    memset(image_array, 0, sizeof(image_array));
    int current_column = half_width - left_pixels;
    int current_row = 0;

    image_channels get_channels;
    // First add the reference to the first row
    FASTA_handler fasta_file(this->ref_file_path);
    string reference_seq = fasta_file.get_reference_sequence(chromosome_name, left_pos, right_pos);
    //reference string
    for(int i=0; i<reference_seq.size(); i++) {
        for(int j=0; j < TOTAL_CHANNELS; j++) {
            if(j==BASE_CHANNEL)
                image_array[current_row][current_column][j] = get_channels.get_base_color(reference_seq[i]);
            else if(j==BASE_QUAL_CHANNEL)
                image_array[current_row][current_column][j] = get_channels.get_base_quality_color(60);
            else if(j==MAP_QUAL_CHANNEL)
                image_array[current_row][current_column][j] = get_channels.get_map_quality_color(40);
            else if(j==STRAND_CHANNEL)
                image_array[current_row][current_column][j] = get_channels.get_strand_color(false);
            else if(j==MATCH_CHANNEL)
                image_array[current_row][current_column][j] = get_channels.get_match_color(true);
            else if(j==CIGAR_CHANNEL)
                image_array[current_row][current_column][j] = get_channels.get_cigar_color(CIGAR_OP_MATCH);
            else if(j==SUPPORT_CHANNEL)
                image_array[current_row][current_column][j] = get_channels.get_support_color(true);
        }
        current_column += 1;
        if(insert_length_map[left_pos+i] > 0){
            for(int l=0; l<insert_length_map[left_pos+i]; l++){
                for(int j=0; j < TOTAL_CHANNELS; j++) {
                    if(j==BASE_CHANNEL)
                        image_array[current_row][current_column][j] = get_channels.get_base_color(reference_seq[i]);
                    else if(j==BASE_QUAL_CHANNEL)
                        image_array[current_row][current_column][j] = get_channels.get_base_quality_color(60);
                    else if(j==MAP_QUAL_CHANNEL)
                        image_array[current_row][current_column][j] = get_channels.get_map_quality_color(40);
                    else if(j==STRAND_CHANNEL)
                        image_array[current_row][current_column][j] = get_channels.get_strand_color(false);
                    else if(j==MATCH_CHANNEL)
                        image_array[current_row][current_column][j] = get_channels.get_match_color(true);
                    else if(j==CIGAR_CHANNEL)
                        image_array[current_row][current_column][j] = get_channels.get_cigar_color(CIGAR_OP_MATCH);
                    else if(j==SUPPORT_CHANNEL)
                        image_array[current_row][current_column][j] = get_channels.get_support_color(true);
                }
                current_column += 1;
            }
        }
    }
    cout<<endl;
    //reference end
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

        // mapping quality
        if(alignment->core.qual <= 0){
            continue;
        }

        // get the base qualities and sequence bases
        uint8_t *seqi = bam_get_seq(alignment);

        // get the cigar operations of the alignment
        uint32_t *cigar = bam_get_cigar(alignment);
        long long read_pos = alignment->core.pos;
        int read_index = 0;
        if(read_pos > right_pos) {
            continue;
        }
        if(read_pos > left_pos) {
            int diff = read_pos - left_pos;
            for(int i=0; i<diff; i++){
                for(int j=0; j < (insert_length_map[left_pos+i] + 1); j++)
                    cout<<" ";
            }
        }

        // initialize a map of read candidates
        map<long long, type_candidate_allele> read_candidate_map;
        map<long long, bool> trace_insert_positions;
        for(int k = 0; k < alignment->core.n_cigar; k++) {
            if(read_pos > right_pos) {
                break;
            }
            int cigar_op = bam_cigar_op(cigar[k]);
            int cigar_len = bam_cigar_oplen(cigar[k]);

            if(cigar_op == BAM_CMATCH){
                if(read_pos + cigar_len < left_pos) {
                    read_pos += cigar_len;
                    read_index += cigar_len;
                    continue;
                }
                //match
                for(int i=0;i<cigar_len;i++) {
                    if(read_pos >= left_pos && read_pos <= right_pos) {
                        if(read_pos - 1 >= alignment->core.pos && insert_length_map[read_pos-1] > 0 && ! trace_insert_positions[read_pos - 1]) {
                            for(int i=0; i<insert_length_map[read_pos-1]; i++) {
                                cout<<"-";
                            }
                        }
                        int ref_index = read_pos - left_pos;
                        cout<<seq_nt16_str[bam_seqi(seqi, read_index)];
                        // get the read base
                        if(seq_nt16_str[bam_seqi(seqi, read_index)] != reference_seq[ref_index]) {
                        }
                    }
                    read_pos += 1;
                    read_index += 1;
                }
            } else if(cigar_op == BAM_CINS) {
                // insert
                long long anchor_position = read_pos - 1;
                trace_insert_positions[anchor_position] = true;
                if(anchor_position >= left_pos && anchor_position < right_pos && read_index > 0) {
                    for (int i = 0; i < cigar_len; i++) {
                        cout<<seq_nt16_str[bam_seqi(seqi, read_index + i)];
                    }
                }
                read_index += cigar_len;
            } else if(cigar_op == BAM_CDEL) {
                // delete
                if(read_pos >= left_pos && read_pos <= right_pos && read_index > 0) {
                    long long delete_len = cigar_len;
                    for(int i=0; i<cigar_len ; i++){
                        if(read_pos >=left_pos && read_pos <= right_pos) {
                            cout<<"*";
                        }
                        read_pos += 1;
                    }
                } else {
                    read_pos += cigar_len;
                }
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
        cout<<endl;
    }

    hts_itr_destroy(iter);
    bam_destroy1(alignment);
}

void image_generator::parse_candidates(string chromosome_name, long long start, long long stop) {
    // find the candidates (labeled or unlabeled)
    candidate_finder candidate_finder_ob;
    map<long long, vector<type_candidate_allele> > positional_candidates;
    map<long long, int> insert_length_map;
    if (this->is_train_mode) {
        positional_candidates = candidate_finder_ob.find_candidates(this->bam_file_path,
                                                                    this->ref_file_path,
                                                                    chromosome_name,
                                                                    start,
                                                                    stop,
                                                                    insert_length_map,
                                                                    this->vcf_file_path,
                                                                    true);
    } else {
        positional_candidates = candidate_finder_ob.find_candidates(this->bam_file_path,
                                                                    this->ref_file_path,
                                                                    chromosome_name,
                                                                    start,
                                                                    stop,
                                                                    insert_length_map);
    }

    for( const auto& pos_info : positional_candidates ) {
        long long pos = pos_info.first;
        vector<type_candidate_allele> candidate_list = pos_info.second;
        for(int i=0;i<candidate_list.size();i++) {
            type_candidate_allele global_candidate = candidate_list[i];
            this->generate_candidate_image(chromosome_name, global_candidate, insert_length_map);
        }
    }
}

image_generator::image_generator(string bam_file_path, string ref_file_path, string vcf_file_path, bool train_mode) {
    this->bam_file_path = bam_file_path;
    this->ref_file_path = ref_file_path;
    this->vcf_file_path = vcf_file_path;
    this->is_train_mode = train_mode;
}

image_generator::~image_generator() {
}