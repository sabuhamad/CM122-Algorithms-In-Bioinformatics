#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <utility>
#include <algorithm>
#include <map>
#include <unordered_map>
#include <set>
#include <thread>
#include <mutex>

using namespace std;
enum class MutationType {
    Substitution,
    Insertion,
    Deletion
};


struct Mutation {
    MutationType type;
    int pos;
    char ref_base;
    char donor_base; //used for substitution & insertion
};


string read_fasta(const string& filename) {
    ifstream input(filename);
    string line, sequence;

    while (getline(input, line)) {
        if (line.empty() || line[0] == '>') {
            continue;
        }
        sequence += line;
    }

    return sequence;
}


vector<pair<string, string>> read_paired_reads(const string& filename) {
    ifstream input(filename);
    if (!input)
        cerr << "ERROR";
    string line;
    vector<pair<string, string>> paired_reads;
    
    size_t estimated_paired_reads = 200000;
    paired_reads.reserve(estimated_paired_reads);
    
    while (getline(input, line)) {
        if (line.empty() || line[0] != '>') {
            continue;
        }
        
        string read1, read2;
        getline(input, read1);
        getline(input, line); //skip second identifier
        getline(input, read2);
        
        paired_reads.push_back({read1, read2});
    }
    
    return paired_reads;
}

unordered_map<string, vector<int>> build_kmer_index(const string& reference, int kmer_size) {
    unordered_map<string, vector<int>> kmer_index;
    for (int i = 0; i < static_cast<int>(reference.size()) - kmer_size + 1; ++i) {
        string kmer = reference.substr(i, kmer_size);
        kmer_index[kmer].push_back(i);
    }
    return kmer_index;
}

int hamming_distance(string read, string reference)
{
    int count = 0;
    for(int i = 0; i < read.length(); i++)
    {
        if(read[i] != reference[i])
        {
            count++;
        }
    }
    return count;
}

void kmer_based_alignment(const string& read, unordered_map<string, vector<int>>& kmer_index, int kmer_size,const  string& reference, vector<unordered_map<char, int>>& base_counts, vector<unordered_map<char, int>>& insertion_map) {
    vector<string> possibleIndel;
    int hamming;
    int align;
    for (int i = 0; i < static_cast<int>(read.size()) - kmer_size + 1; i += kmer_size) {
        bool detected = false;
        bool perfect = false;
        
        string kmer = read.substr(i, kmer_size);
        if(kmer_index[kmer].empty())
        {
            continue;
        }
        if(kmer_index[kmer][0] < i)
        {
            continue;
        }
        hamming = hamming_distance(reference.substr(kmer_index[kmer][0], read.length()), read);
        align = kmer_index[kmer][0] - i;
        for(auto pos : kmer_index[kmer])
        {
            int distance = hamming_distance(reference.substr(pos - i, read.length()), read);
            if(distance < hamming)
            {
                hamming = distance;
                align = pos - i;
            }
        }
        
        if(hamming < 3)
        {
            detected = true;
            if(read == reference.substr(align, read.length()))
            {
                perfect = true;
                break;
            }
            
            for(int j = 0; j < read.length(); j++)
            {
                if(read[j] != reference[align + j])
                {
                    base_counts[align+j][read[j]]++;
                }
            }
        }
        if(detected != true)
        {
            possibleIndel.push_back(read);
        }
    }
    
    for (auto indel : possibleIndel)
    {
        for (int i = 0; i < indel.length()- kmer_size+1; i += kmer_size){
            string kmer = indel.substr(i,kmer_size);

            if (kmer_index[kmer].empty())
                continue;
            for (auto pos : kmer_index[kmer]){
                bool detected = false;
                if (pos-i < 0 || pos-i+indel.length() >= reference.length())
                    continue;
                string alignment = reference.substr(pos-i, indel.length());
                for (int b = 0; b < indel.length(); b++){
                    if (indel[b] != alignment[b]){
                        if (pos-i+b+1 < 0 || pos-i+indel.length()+1 >= reference.length())
                            continue;
                        string substring_of_read = indel.substr(b, indel.length()-b);
                        string align_substring = reference.substr(pos-i+b+1, indel.length()-b);
                        if (substring_of_read == align_substring){
                            base_counts[pos - i + b - 1]['-']++;
                            detected = true;
                            break;
                        }
                        substring_of_read = indel.substr(b+1, indel.length()-b-1);
                        align_substring = reference.substr(pos-i+b, indel.length()-b-1);
                        if (substring_of_read == align_substring){
                            base_counts[pos-i+b - 1]['+']++;
                            insertion_map[pos-i+b - 1][indel[b]]++;
                            detected = true;
                            break;
                        }
                    }
                }
                if (detected == true)
                {
                    break;
                }
            }
        }
    }
    
}

vector<Mutation> identify_mutations(const string& reference, const vector<pair<string, string>>& paired_reads, double mutation_threshold, int kmer_size) {
    cout << "Identifying mutations..." << endl;

    vector<unordered_map<char, int>> base_counts(reference.size());
    vector<unordered_map<char, int>> insertion_map(reference.size());


    unordered_map<string, vector<int>> kmer_index = build_kmer_index(reference, kmer_size);

    for (const auto& read_pair : paired_reads) {
        for (const auto& read : {read_pair.first, read_pair.second}) {
            if(read.size() >= kmer_size)
            {
                kmer_based_alignment(read, kmer_index, kmer_size, reference, base_counts, insertion_map);
            }
        }
    }

    vector<Mutation> mutations;

    for (int pos = 0; pos < reference.size(); ++pos) {
        int total_count = 0;
        for (const auto& base_count : base_counts[pos]) {
            total_count += base_count.second;
        }

        char ref_base = reference[pos];
        for (const auto& base_count : base_counts[pos]) {
            char base = base_count.first;
            int count = base_count.second;

            double mutation_frequency = (double)count / (double)(total_count);
            if (base != ref_base && (mutation_frequency >= mutation_threshold && count >= 3)) {
                Mutation mutation;

                if (base == '+') {
                    mutation.type = MutationType::Insertion;
                } else if (base == '-') {
                    mutation.type = MutationType::Deletion;
                } else {
                    mutation.type = MutationType::Substitution;
                }

                mutation.pos = pos;
                mutation.ref_base = ref_base;

                if (base == '+') {
                    
                    char max_inserted_base = ' ';
                    int max_inserted_count = 0;
                    for (const auto& donor_base : insertion_map[pos]) {
                        if (donor_base.second > max_inserted_count) {
                            max_inserted_base = donor_base.first;
                            max_inserted_count = donor_base.second;
                        }
                    }
                    mutation.donor_base = max_inserted_base;
                } else if (base == '-') {
                    mutation.donor_base = ' ';
                } else {
                    mutation.donor_base = base;
                }

                mutations.push_back(mutation);
            }
        }
    }

    return mutations;
}


int main() {
    string ref_fasta = "project1b_1000000_reference_genome.fasta";
    string ref_paired_reads = "project1b_1000000_with_error_paired_reads.fasta";
    string reference = read_fasta(ref_fasta);
    vector<pair<string, string>> paired_reads = read_paired_reads(ref_paired_reads);

    double mutation_threshold = 0.5;
    int kmer_size = 15;
    //identify mutations
    vector<Mutation> mutations = identify_mutations(reference, paired_reads, mutation_threshold, kmer_size);

    //print mutations
    for (const auto& mutation : mutations) {
        
        if (mutation.type == MutationType::Substitution) {
            cout << ">S" << mutation.pos << " ";
            cout << mutation.ref_base << " " << mutation.donor_base << endl;
        }
    }
    for(const auto& mutation : mutations)
    {
        if(mutation.type == MutationType::Insertion){
            cout << ">I" << mutation.pos << " ";
            cout << mutation.donor_base << endl;
        }
    }
    for(const auto& mutation : mutations)
    {
        if(mutation.type == MutationType::Deletion)
        {
            cout << ">D" << mutation.pos << " ";
            cout << mutation.ref_base;
            cout << endl;
        }
    }

    return 0;
}
