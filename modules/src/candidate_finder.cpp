//
// Created by Kishwar Shafin on 6/18/18.
//
#include "../headers/candidate_finder.h"

int candidate_finder::get_genotype(vector<int> genotype) {
    if(genotype.size() == 0) {
        return -1;
    }
    else if(genotype.size() < 2) {
        return genotype[0];
    }
    else if(genotype.size() == 2) {
        if(genotype[0] == 0 && genotype[1] == 0) return GENOTYPE_HOM;
        else if(genotype[0] == 0 && genotype[1] == 1) return GENOTYPE_HET;
        else if(genotype[0] == 1 && genotype[1] == 0) return GENOTYPE_HET;
        else if(genotype[0] == 0 && genotype[1] == 2) return GENOTYPE_HET;
        else if(genotype[0] == 2 && genotype[1] == 0) return GENOTYPE_HET;
        else if(genotype[0] == 1 && genotype[1] == 1) return GENOTYPE_HOMALT;
        else if(genotype[0] == 2 && genotype[1] == 2) return GENOTYPE_HOMALT;
        else if(genotype[0] == 1 && genotype[1] == 2) return GENOTYPE_HET;
        else if(genotype[0] == 2 && genotype[1] == 1) return GENOTYPE_HET;
        return -2;
    }
    return -3;
}


map<long long, vector<type_candidate_allele> > candidate_finder::find_candidates(string bam_file_path,
                                                                                 string ref_file_path,
                                                                                 string chromosome,
                                                                                 long long start,
                                                                                 long long stop,
                                                                                 map<long long, int> &insert_length_map,
                                                                                 string vcf_file_path,
                                                                                 bool label_candidates) {
    BAM_handler bam_file(bam_file_path);
    FASTA_handler fasta_file(ref_file_path);
    string reference_seq = fasta_file.get_reference_sequence(chromosome, start, stop);

    // get the id of the chromosome or sequence name
    const int tid = bam_name2id(bam_file.header, chromosome.c_str());

    // get the iterator
    hts_itr_t *iter  = sam_itr_queryi(bam_file.idx, tid, start, stop);

    // initialize an alignment
    bam1_t* alignment = bam_init1();
    map<type_candidate_allele, int> global_candidate_map;

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

        // initialize a map of read candidates
        map<long long, vector<type_candidate_allele> > read_candidate_map;

        for(int k = 0; k < alignment->core.n_cigar; k++) {
            int cigar_op = bam_cigar_op(cigar[k]);
            int cigar_len = bam_cigar_oplen(cigar[k]);

            if(cigar_op == BAM_CMATCH){
                //match
                for(int i=0;i<cigar_len;i++) {
                    if(read_pos >= start && read_pos <= stop) {
                        int ref_index = read_pos - start;
                        // get the read base
                        if(seq_nt16_str[bam_seqi(seqi, read_index)] != reference_seq[ref_index]) {
                            type_candidate_allele candidate;
                            candidate.pos = read_pos;
                            candidate.allele = seq_nt16_str[bam_seqi(seqi, read_index)];
                            candidate.candidate_type = SNP_TYPE;
                            read_candidate_map[read_pos].push_back(candidate);
                        }
                    }
                    read_pos += 1;
                    read_index += 1;
                }
            } else if(cigar_op == BAM_CINS) {
                // insert
                long long anchor_position = read_pos - 1;
                if(anchor_position >= start && anchor_position <= stop && read_index > 0) {
                    string insert_allele;
                    for (int i = -1; i < cigar_len; i++) {
                        insert_allele += seq_nt16_str[bam_seqi(seqi, read_index + i)];
                    }
                    type_candidate_allele candidate;
                    candidate.pos = anchor_position;
                    candidate.allele = insert_allele;
                    candidate.candidate_type = INSERT_TYPE;
                    read_candidate_map[anchor_position].push_back(candidate);
                }
                if(read_index > 0) {
                    if(insert_length_map.find(anchor_position) != insert_length_map.end())
                        insert_length_map[anchor_position] = max(insert_length_map[anchor_position], cigar_len);
                    else
                        insert_length_map[anchor_position] = cigar_len;
                }
                read_index += cigar_len;
            } else if(cigar_op == BAM_CDEL) {
                // delete
                long long anchor_position = read_pos - 1;
                if(anchor_position >= start && anchor_position + cigar_len <= stop && read_index > 0) {
                    // get the anchor allele from the read
                    string delete_allele;
                    delete_allele = seq_nt16_str[bam_seqi(seqi, read_index-1)];
                    int ref_index = read_pos - start;
                    for (int i = 0; i < cigar_len; i++) {
                        delete_allele += reference_seq[ref_index + i];
                    }
                    type_candidate_allele candidate;
                    candidate.pos = anchor_position;
                    candidate.allele = delete_allele;
                    candidate.candidate_type = DELETE_TYPE;
                    // cout<<"Delete: "<<read_pos<<" "<<delete_allele<<endl;
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

        for( const auto& read_candidate_it : read_candidate_map ) {
            vector<type_candidate_allele> candidates = read_candidate_it.second;
            for(int i=0; i < candidates.size(); i++) {
                global_candidate_map[candidates[i]] += 1;
            }
        }
    }

    map<type_candidate_allele, bool> positional_top_candidates;
    for( const auto& candidate_info : global_candidate_map ) {
        type_candidate_allele candidate = candidate_info.first;
        int freq = candidate_info.second;
        candidate.freq = freq;
        positional_top_candidates[candidate] = true;
    }

    map<long long, vector<type_vcf_record> > positional_variants;
    if(label_candidates){
        vcf_handler vcf_file(vcf_file_path);
        positional_variants = vcf_file.get_positional_vcf_records(chromosome, start - 50, stop + 50);
    }
    map<long long, vector<type_candidate_allele> > positional_candidates;

    for( const auto& candidate_info : positional_top_candidates ) {
        type_candidate_allele candidate = candidate_info.first;

        if(candidate.freq <= 1 && candidate.candidate_type==SNP_TYPE) {
            continue;
        }
        if(positional_candidates[candidate.pos].size() < 2){
            int candidate_gt = GENOTYPE_HOM;
            bool is_labeled = false;
            if(label_candidates) {
                is_labeled = true;
                if(positional_variants[candidate.pos].size() > 0){
                    vector<type_vcf_record> vcf_records = positional_variants[candidate.pos];
                    // for each record in that position
                    for(int i=0; i<vcf_records.size(); i++) {
                        // if there's an alt allele
                        for(int j=0;j<vcf_records[i].alt_allele.size();j++) {
                            // that matches the candidate
                            if(vcf_records[i].alt_allele[j].alt_type == candidate.candidate_type &&
                               vcf_records[i].alt_allele[j].alt_allele == candidate.allele &&
                               vcf_records[i].is_filter_pass){
                                // then find the genotype
                                candidate_gt = this->get_genotype(vcf_records[i].genotype);
                                if(candidate_gt < 0) {
                                    if(candidate_gt == -1){
                                        cerr<<"WARN: EMPTY GENOTYPE VECTOR."<<vcf_records[i].chromosome_name<<" ";
                                        cerr<<vcf_records[i].start_pos<<" "<<vcf_records[i].id<<endl;
                                    } else {
                                        cerr<<"WARN: INVALID GENOTYPE FOUND."<<vcf_records[i].chromosome_name<<" ";
                                        cerr<<vcf_records[i].start_pos<<" "<<vcf_records[i].id<<endl;
                                    }
                                }
                            }
                        }
                    }
                }
            }
            candidate.is_labeled = is_labeled;
            candidate.labeled_genotype = candidate_gt;
            positional_candidates[candidate.pos].push_back(candidate);
        }
    }

    hts_itr_destroy(iter);
    bam_destroy1(alignment);

    return positional_candidates;
}

candidate_finder::candidate_finder() {
}

candidate_finder::~candidate_finder() {
}