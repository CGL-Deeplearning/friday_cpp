#include <stdio.h>
#include <iostream>
#include <string.h>
#include <vector>
#include "modules/headers/bam_handler.h"
#include "modules/headers/fasta_handler.h"
#include "modules/headers/vcf_handler.h"
#include "modules/headers/candidate_finder.h"
using namespace std;

int main(int argc, char **argv)
{
    string bam_file_path = "/data/users/common/GIAB/bam/chr19_GRCh37.bam";
    string fasta_file_path = "/data/users/common/GIAB/ref/GRCh37_WG.fa";
    string vcf_file_path = "/data/users/common/GIAB/vcf/GRCh37_vcf/GRCh37_GIAB.vcf.gz";
    string region = "19";
    int  start_pos = 5356136;
    int end_pos = 7275660;


    // Sample names
    /*set<string> sample_names;
    sample_names = bam_file.get_sample_names();

    // chromosome names
    vector <type_sequence> chromosome_names;
    chromosome_names = bam_file.get_chromosome_sequence_names();*/

    // get reference seq


    candidate_finder candidate_finder_oj;
    map<long long, vector<type_candidate_allele> > positional_candidates;
    positional_candidates = candidate_finder_oj.find_candidates(bam_file_path,
                                                                fasta_file_path,
                                                                region,
                                                                start_pos,
                                                                end_pos,
                                                                vcf_file_path,
                                                                false);

    for( const auto& pos_info : positional_candidates ) {
        long long pos = pos_info.first;
        vector<type_candidate_allele> candidate_list = pos_info.second;
        for(int i=0;i<candidate_list.size();i++) {
            type_candidate_allele global_candidate = candidate_list[i];
            cout << pos << " " << global_candidate.pos << " " << global_candidate.allele << " ";
            cout << global_candidate.candidate_type << " " << global_candidate.freq << " ";
            cout << global_candidate.labeled_genotype << endl;
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
