#include <stdio.h>
#include <iostream>
#include <string.h>
#include <vector>
#include "modules/headers/train_data_generator.h"
#include "modules/headers/bed_handler.h"
#include "tqdm.h"
using namespace std;

int main(int argc, char **argv)
{
    string bam_file_path = "/data/users/common/GIAB/bam/chr19_GRCh37.bam";
    string ref_file_path = "/data/users/common/GIAB/ref/GRCh37_WG.fa";
    string vcf_file_path = "/data/users/common/GIAB/vcf/GRCh37_vcf/GRCh37_GIAB.vcf.gz";
    string confident_bed_path = "/data/users/common/GIAB/confident_bed/GRCh37_confident/NA12878_GRCh37_confident.bed";


    bool is_JARVIS = 1;

    if(is_JARVIS) {
        bam_file_path = "/data/users/common/GIAB/bam/NA12878_GIAB_30x_GRCh37.sorted.bam";
        ref_file_path = "/data/users/common/GIAB/ref/GRCh37_WG.fa";
        vcf_file_path = "/data/users/common/GIAB/vcf/GRCh37_vcf/NA12878_GRCh37.vcf.gz";
        confident_bed_path = "/data/users/common/GIAB/confident_bed/GRCh37_confident/NA12878_GRCh37_confident.bed";
    }


    train_data_generator image_generator(bam_file_path,
                                         ref_file_path,
                                         vcf_file_path,
                                         confident_bed_path);
    image_generator.genome_level_processes();
    // Sample names
    /*set<string> sample_names;
    sample_names = bam_file.get_sample_names();

    // chromosome names
    vector <type_sequence> chromosome_names;
    chromosome_names = bam_file.get_chromosome_sequence_names();*/


    /*for(int i=0; i<candidate_positions.size(); i++) {
        long long pos = candidate_positions[i];
        vector<type_candidate_allele> candidate_list = candidate_lists[i];

        vector<image_generator> generated_images;
        for(int i=0; i< candidate_list.size(); i++) {
            // candidate allele i
            image_generator image_generator_ob(chromosome_name, bam_file_path, ref_file_path,
                                               candidate_list[i], insert_length_map);
            image_generator_ob.generate_candidate_image();

            generated_images.push_back(image_generator_ob);

            string output_file_name = output_directory + chromosome_name + "_" +
                                      to_string(candidate_list[i].pos) + "_" + to_string(i) + ".h5";
            hdf5_handler image_saver(output_file_name);
            image_saver.save_img_to_hdf5(image_generator_ob.image_array);
            image_saver.save_allele_info(chromosome_name, to_string(candidate_list[i].pos),
                                         candidate_list[i].allele, to_string(candidate_list[i].candidate_type),
                                         ".", "-1");
        }
        if(candidate_list.size() == 2) {
            image_generator ig;
            ig.generate_combined_image(generated_images[0], generated_images[1]);

            string output_file_name = output_directory + chromosome_name + "_" +
                                      to_string(candidate_list[0].pos) + "_2" + ".h5";

            hdf5_handler image_saver(output_file_name);
            image_saver.save_img_to_hdf5(ig.image_array);
            image_saver.save_allele_info(chromosome_name, to_string(candidate_list[0].pos),
                                         candidate_list[0].allele, to_string(candidate_list[0].candidate_type),
                                         candidate_list[1].allele, to_string(candidate_list[1].candidate_type));

        }
        progress_bar.progress(progress, total_candidates);
        progress += 1;
    }*/
    //cout<<reference_seq<<endl;

    /*for(int i=0; i<reads.size(); i++){
        string read_seqence = reads[i].sequence;
        int read_start_pos = reads[i].pos;
        string reference_seq = FASTA_file.get_reference_sequence(region, reads[i].pos, reads[i].pos+read_seqence.size()-1);
        //cout<<"READ:\t"<<read_seqence<<endl;
        //cout<<"REFS:\t"<<reference_seq<<endl;
    }*/

    /*string vcf_file_path = "/data/users/common/GIAB/vcf/GRCh37_vcf/GRCh37_GIAB.vcf.gz";
    vcf_handler vcf_file(vcf_file_path);
    vector<type_vcf_record> vcf_records = vcf_file.get_vcf_records(region, start_pos, end_pos);

    cout<<"Chr\tSt\tend\tqual\tphased\tsn\tgt\tref\t\talts\t\ttype"<<endl;
    for(int i=0;i<vcf_records.size(); i++){
        cout<<vcf_records[i].chromosome_name<<"\t"<<vcf_records[i].start_pos<<"\t"<<vcf_records[i].end_pos<<"\t";
        cout<<vcf_records[i].qual<<"\t"<<vcf_records[i].is_phased<<"\t"<<vcf_records[i].sample_name<<"\t";
        for(int j=0;j<vcf_records[i].genotype.size();j++){
            cout<<vcf_records[i].genotype[j]<<",";
        }
        cout<<"\t";
        for(int j=0;j<vcf_records[i].alt_allele.size();j++){
            cout<<vcf_records[i].alt_allele[j].ref<<",";
        }
        cout<<"\t";
        for(int j=0;j<vcf_records[i].alt_allele.size();j++){
            cout<<vcf_records[i].alt_allele[j].alt_allele<<",";
        }
        cout<<"\t";
        for(int j=0;j<vcf_records[i].alt_allele.size();j++){
            cout<<vcf_records[i].alt_allele[j].alt_type<<",";
        }
        cout<<endl;
    }*/
	return 0;
}
