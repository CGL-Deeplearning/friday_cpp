#include <stdio.h>
#include <iostream>
#include <string.h>
#include <vector>
#include "modules/headers/candidate_finder.h"
#include "modules/headers/image_generator.h"
#include "modules/headers/hdf5_handler.h"
using namespace std;

int main(int argc, char **argv)
{
    string bam_file_path = "/data/users/common/GIAB/bam/chr19_GRCh37.bam";
    string ref_file_path = "/data/users/common/GIAB/ref/GRCh37_WG.fa";
    string vcf_file_path = "/data/users/common/GIAB/vcf/GRCh37_vcf/GRCh37_GIAB.vcf.gz";
    string chromosome_name = "19";
    long long  start = 5356142;
    long long stop = 5356148;
    bool is_train_mode = true;


    // Sample names
    /*set<string> sample_names;
    sample_names = bam_file.get_sample_names();

    // chromosome names
    vector <type_sequence> chromosome_names;
    chromosome_names = bam_file.get_chromosome_sequence_names();*/

    // find the candidates (labeled or unlabeled)
    candidate_finder candidate_finder_ob;
    map<long long, vector<type_candidate_allele> > positional_candidates;
    map<long long, int> insert_length_map;

    positional_candidates = candidate_finder_ob.find_candidates(bam_file_path,
                                                                ref_file_path,
                                                                chromosome_name,
                                                                start,
                                                                stop,
                                                                insert_length_map,
                                                                vcf_file_path,
                                                                true);

    for( const auto& pos_info : positional_candidates ) {
        long long pos = pos_info.first;
        vector<type_candidate_allele> candidate_list = pos_info.second;
        for(int i=0;i<candidate_list.size();i++) {
            cout<<candidate_list[i].pos<<" "<<candidate_list[i].allele<<" "<<candidate_list[i].candidate_type<<endl;
            image_generator image_generator_ob(chromosome_name, bam_file_path, ref_file_path, candidate_list[i], insert_length_map);
            image_generator_ob.generate_candidate_image();
            hdf5_handler image_saver;
            image_saver.save_img_to_hdf5(image_generator_ob);
            exit(0);
        }
    }
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
