#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <utility>
#include <algorithm>
#include <map>
#include <unordered_map>
#include <set>

using namespace std;

enum class MutationType {
    Substitution,
    Insertion,
    Deletion
};

struct Mutation {
    MutationType type;
    int position;
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

//implement smith waterman alignment
pair<int, vector<Mutation>> smith_waterman_alignment(const string& read, const string& reference, int match_score, int mismatch_penalty, int open_penalty, int up_score) {
    int n = read.size();
    int m = reference.size();
    vector<vector<int>> dp(n + 1, vector<int>(m + 1, 0));
    vector<vector<pair<int, int>>> traceback(n + 1, vector<pair<int, int>>(m + 1));

    int max_score = 0;
    pair<int, int> max_score_position;

    for (int i = 1; i <= n; ++i) {
        for (int j = 1; j <= m; ++j) {
            int match_mismatch_score = dp[i - 1][j - 1] + (read[i - 1] == reference[j - 1] ? match_score : mismatch_penalty);
            int deletion_score = max(dp[i - 1][j] + open_penalty, dp[i - 1][j] + up_score);
            int insertion_score = max(dp[i][j - 1] + open_penalty, dp[i][j - 1] + up_score);

            dp[i][j] = max({0, match_mismatch_score, deletion_score, insertion_score});
            if (dp[i][j] == match_mismatch_score) {
                traceback[i][j] = {i - 1, j - 1};
            } else if (dp[i][j] == deletion_score) {
                traceback[i][j] = {i - 1, j};
            } else {
                traceback[i][j] = {i, j - 1};
            }

            if (dp[i][j] > max_score) {
                max_score = dp[i][j];
                max_score_position = {i, j};
            }
        }
    }
    //store all possible mutations
    vector<Mutation> mutations;
    int i = max_score_position.first;
    int j = max_score_position.second;
    while (i > 0 && j > 0) {
        int prev_i = traceback[i][j].first;
        int prev_j = traceback[i][j].second;

        if (prev_i == i - 1 && prev_j == j - 1) {
            if (read[i - 1] != reference[j - 1]) {
                mutations.push_back({MutationType::Substitution, j - 1, reference[j - 1], read[i - 1]});
            }
        } else if (prev_i == i - 1) {
            mutations.push_back({MutationType::Insertion, j - 1, '+', read[i - 1]});
        } else {
            mutations.push_back({MutationType::Deletion, j - 1, reference[j - 1], '-'});
        }



        i = prev_i;
        j = prev_j;
    }

    reverse(mutations.begin(), mutations.end());

    return {max_score, mutations};
}



vector<Mutation> identify_mutations(const string& reference, const vector<pair<string, string>>& paired_reads, double mutation_threshold) {
    vector<Mutation> mutations;

    //alignment parameters
    int match_score = 2;
    int mismatch_penalty = -1;
    int open_penalty = -2;
    int up_score = -2;

    vector<unordered_map<char, int>> base_counts(reference.size());
    vector<unordered_map<char,int>> insertion_map(reference.size());
    for (const auto& read_pair : paired_reads) {
        for (const auto& read : {read_pair.first, read_pair.second}) {
            pair<int, vector<Mutation>> alignment_result = smith_waterman_alignment(read, reference, match_score, mismatch_penalty, open_penalty, up_score);
            vector<Mutation> read_mutations = alignment_result.second;

            for (const auto& mutation : read_mutations) {
                if (mutation.type == MutationType::Insertion) {
                    base_counts[mutation.position]['+']++;
                    insertion_map[mutation.position][mutation.donor_base]++;
                } else if (mutation.type == MutationType::Deletion) {
                    base_counts[mutation.position]['-']++;
                } else {
                    base_counts[mutation.position][mutation.donor_base]++;
                }
            }
        }
    }

    // Identify mutations
 
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
                mutation.position = pos;
                mutation.ref_base = ref_base;
                if (base == '+')
                {
                    for(const auto& donor_base : insertion_map[pos])
                        mutation.donor_base = donor_base.first;
                }
                else if(base == '-') {
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
    string ref_fasta = "/Users/sabuhamad/Desktop/CM122/Project1/sample_1000/sample_1000_reference_genome.fasta";
    string ref_paired_reads = "/Users/sabuhamad/Desktop/CM122/Project1/sample_1000/sample_1000_with_error_single_reads.fasta";
    string reference = read_fasta(ref_fasta);
    vector<pair<string, string>> paired_reads = read_paired_reads(ref_paired_reads);

    double mutation_threshold = 0.4;
    
    //identify mutations
    vector<Mutation> mutations = identify_mutations(reference, paired_reads, mutation_threshold);

    //print mutations
    for (const auto& mutation : mutations) {
        
        if (mutation.type == MutationType::Substitution) {
            cout << ">S" << mutation.position << " ";
            cout << mutation.ref_base << " " << mutation.donor_base << endl;
        }
    }
    for(const auto& mutation : mutations)
    {
        if(mutation.type == MutationType::Insertion){
            cout << ">I" << mutation.position << " ";
            cout << mutation.donor_base << endl;
        }
    }
    for(const auto& mutation : mutations)
    {
        if(mutation.type == MutationType::Deletion)
        {
            cout << ">D" << mutation.position << " ";
            cout << mutation.ref_base;
            cout << endl;
        }
    }

    return 0;
}
