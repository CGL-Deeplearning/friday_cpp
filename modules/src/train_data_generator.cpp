//
// Created by Kishwar Shafin on 7/5/18.
//
#include "../headers/train_data_generator.h"
train_data_generator::train_data_generator(string bam_file_path,
                                           string ref_file_path,
                                           string vcf_file_path,
                                           string confident_bed_path) {
    this->chromosome_name = chromosome_name;
    this->bam_file_path = bam_file_path;
    this->ref_file_path = ref_file_path;
    this->vcf_file_path = vcf_file_path;
    this->confident_bed_path = confident_bed_path;
}

void train_data_generator::generate_labeled_images(string chromosome_name, long long start_pos, long long end_pos) {
    // find the candidates (labeled or unlabeled)
    candidate_finder candidate_finder_ob;
    map<long long, vector<type_candidate_allele> > positional_candidates;
    map<long long, int> insert_length_map;

    positional_candidates = candidate_finder_ob.find_candidates(bam_file_path,
                                                                ref_file_path,
                                                                chromosome_name,
                                                                start_pos,
                                                                end_pos,
                                                                insert_length_map,
                                                                vcf_file_path,
                                                                true);

    vector<long long> candidate_positions;
    vector<vector <type_candidate_allele> > candidate_lists;

    for( const auto& pos_info : positional_candidates ) {
        if (pos_info.second.size() == 0) continue;
        candidate_positions.push_back(pos_info.first);
        candidate_lists.push_back(pos_info.second);
    }

    string output_directory = "./outputs/";

//    tqdm progress_bar;
//    int progress = 0;
//    int total_ite = candidate_positions.size();

    #pragma omp parallel for
    for(int j=0; j<candidate_positions.size(); j++) {
        vector<image_generator> generated_images;

        #pragma omp parallel for
        for(int i=0; i< candidate_lists[j].size(); i++) {
            // candidate allele i
            image_generator image_generator_ob(chromosome_name, bam_file_path, ref_file_path,
                                               candidate_lists[j][i], insert_length_map);
            image_generator_ob.generate_candidate_image();

            generated_images.push_back(image_generator_ob);

            string output_file_name = output_directory + chromosome_name + "_" +
                                      to_string(candidate_lists[j][i].pos) + "_" + to_string(i) + ".h5";
            #pragma omp critical
            {
                hdf5_handler image_saver(output_file_name);
                image_saver.save_img_to_hdf5(image_generator_ob.image_array);
                image_saver.save_allele_info(chromosome_name, to_string(candidate_lists[j][i].pos),
                                             candidate_lists[j][i].allele, to_string(candidate_lists[j][i].candidate_type),
                                             ".", "-1");
            }
        }
        if(candidate_lists[j].size() == 2) {
            image_generator ig;
            ig.generate_combined_image(generated_images[0], generated_images[1]);
            string output_file_name = output_directory + chromosome_name + "_" +
                                      to_string(candidate_lists[j][0].pos) + "_2" + ".h5";
            #pragma omp critical
            {
                hdf5_handler image_saver(output_file_name);
                image_saver.save_img_to_hdf5(ig.image_array);
                image_saver.save_allele_info(chromosome_name, to_string(candidate_lists[j][0].pos),
                                             candidate_lists[j][0].allele, to_string(candidate_lists[j][0].candidate_type),
                                             candidate_lists[j][1].allele, to_string(candidate_lists[j][1].candidate_type));
            }
        }
//        progress += 1;
//        progress_bar.progress(progress, total_ite);
    }
}

void train_data_generator::genome_level_processes() {
    bed_handler bed_processor(confident_bed_path);
    bed_processor.read_bed_file();

    BAM_handler bam_processor(bam_file_path);

    // get all the sequence names
    vector <string> bam_names;
    bam_names = bam_processor.get_chromosome_sequence_names();

    vector<string> filtered_sequence_names;
    for(int i=0; i < bed_processor.chromosome_name_set.size(); i++) {
        string bed_chromosome_name = bed_processor.chromosome_name_set[i];
        if(find(bam_names.begin(), bam_names.end(), bed_chromosome_name) != bam_names.end())
            filtered_sequence_names.push_back(bed_chromosome_name);
        else
            printf(ANSI_COLOR_YELLOW "WARN: BED SEQUENCE NOT FOUND IN BAM FILE %s\n" ANSI_COLOR_RESET,
                   bed_chromosome_name.c_str());
    }
    printf(ANSI_COLOR_BLUE "INFO: TOTAL SEQUENCES FOUND %d\n" ANSI_COLOR_RESET, int(filtered_sequence_names.size()));

    vector<type_sequence> sequence_with_length;
    sequence_with_length = bam_processor.get_chromosome_sequence_names_with_length();
    map<string, long long> length_by_sequence;
    for(int i=0; i < sequence_with_length.size(); i++) {
        type_sequence seq = sequence_with_length[i];
        length_by_sequence[seq.sequence_name] = seq.sequence_length;
    }

    vector<string> sequences = {"19"};
    unsigned int available_threads = thread::hardware_concurrency();

    long long each_segment_length;

    for(int i=0; i<sequences.size(); i++) {
        string chromosome_name = sequences[i];
        long long len = length_by_sequence[chromosome_name];
        long long current_pos = 0;
        each_segment_length = len / available_threads;
        vector<bed_interval> intervals;
        while(current_pos < len) {
            bed_interval in(current_pos, min(len, current_pos+each_segment_length));
            current_pos += each_segment_length;

            intervals.push_back(in);
        }

        printf(ANSI_COLOR_GREEN "INFO: GENERATING IMAGES FROM %s\n" ANSI_COLOR_RESET, chromosome_name.c_str());
        tqdm progress_bar;
        int progress = 0;
        int total_ite = intervals.size();

        progress_bar.set_label(chromosome_name);
        #pragma omp parallel
        {
            for (int j = 0; j < total_ite; j++) {
                progress_bar.progress(progress, total_ite);
                progress += 1;
                generate_labeled_images(chromosome_name, intervals[j].start_pos, intervals[j].end_pos);
            }
        }

        printf(ANSI_COLOR_BLUE "\nINFO: IMAGE GENERATION FINISHED %s\n" ANSI_COLOR_RESET, chromosome_name.c_str());
    }
}

train_data_generator::~train_data_generator() {
}