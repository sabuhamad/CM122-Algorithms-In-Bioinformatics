#include <fstream>
#include <iostream>
#include <string>
#include <unordered_map>
#include <sstream>
#include <vector>

using namespace std;

vector<vector<double>> read_pwm_file(const string &filename) {
  vector<vector<double>> pwm;
  ifstream file(filename);

  if (file.is_open()) {
    string line;

    while (getline(file, line)) {
      stringstream ss(line);
      double value;
      vector<double> row;

      while (ss >> value) {
        row.push_back(value);
        if (ss.peek() == '\t') {
          ss.ignore();
        }
      }

      pwm.push_back(row);
    }

    file.close();
  }

  return pwm;
}


vector<vector<double>> generate_pwm_from_motifs(const vector<string> &motifs, double pseudo_count = 1) {
  vector<vector<double>> pwm(4, vector<double>(motifs[0].size(), pseudo_count));

  for (const auto &motif : motifs) {
    for (size_t i = 0; i < motif.size(); ++i) {
      switch (motif[i]) {
        case 'A': ++pwm[0][i]; break;
        case 'C': ++pwm[1][i]; break;
        case 'G': ++pwm[2][i]; break;
        case 'T': ++pwm[3][i]; break;
      }
    }
  }

  for (size_t i = 0; i < pwm[0].size(); ++i) {
    double col_sum = pwm[0][i] + pwm[1][i] + pwm[2][i] + pwm[3][i];
    pwm[0][i] /= col_sum;
    pwm[1][i] /= col_sum;
    pwm[2][i] /= col_sum;
    pwm[3][i] /= col_sum;
  }

  return pwm;
}


double prob_kmer_given_pwm(const string &kmer, const vector<vector<double>> &pwm) {
  double prob = 1.0;
  for (size_t i = 0; i < kmer.size(); ++i) {
    switch (kmer[i]) {
      case 'A': prob *= pwm[0][i]; break;
      case 'C': prob *= pwm[1][i]; break;
      case 'G': prob *= pwm[2][i]; break;
      case 'T': prob *= pwm[3][i]; break;
    }
  }
  return prob;
}


double calculate_pwm_score(const string &sequence, const vector<vector<double>> &pwm) {
  double max_score = -numeric_limits<double>::infinity();
  int motif_length = pwm[0].size();

  for (size_t i = 0; i <= sequence.size() - motif_length; ++i) {
    double score = 0;

    for (int j = 0; j < motif_length; ++j) {
      char nucleotide = sequence[i + j];

      switch (nucleotide) {
        case 'A':
          score += pwm[0][j];
          break;
        case 'C':
          score += pwm[1][j];
          break;
        case 'G':
          score += pwm[2][j];
          break;
        case 'T':
          score += pwm[3][j];
          break;
        default:
          break;
      }
    }

    max_score = max(max_score, score);
  }

  return max_score;
}

vector<vector<double>> greedy_motif_search(const vector<string> &sequences, int k, int num_iterations = 3) {

  vector<string> best_motifs;
  double best_pwm_score = -numeric_limits<double>::infinity();
  vector<vector<double>> best_pwm;

  for (int iter = 0; iter < num_iterations; ++iter) {
    vector<string> motifs;
    for (const auto &seq : sequences) {
      int start_idx = rand() % (seq.size() - k + 1);
      motifs.push_back(seq.substr(start_idx, k));
    }

    bool converged = false;
    while (!converged) {
      converged = true;
      vector<vector<double>> pwm = generate_pwm_from_motifs(motifs);

      for (size_t i = 0; i < sequences.size(); ++i) {
        double max_prob = -1;
        string best_kmer;

        for (size_t j = 0; j <= sequences[i].size() - k; ++j) {
          string kmer = sequences[i].substr(j, k);
          double prob = prob_kmer_given_pwm(kmer, pwm);

          if (prob > max_prob) {
            max_prob = prob;
            best_kmer = kmer;
          }
        }

        if (motifs[i] != best_kmer) {
          motifs[i] = best_kmer;
          converged = false;
        }
      }
    }

    double pwm_score = 0;
    vector<vector<double>> pwm = generate_pwm_from_motifs(motifs);
    for (const auto &motif : motifs) {
      pwm_score += calculate_pwm_score(motif, pwm);
    }

    if (pwm_score > best_pwm_score) {
      best_pwm_score = pwm_score;
      best_motifs = motifs;
      best_pwm = pwm;
    }
  }

  return best_pwm;
}


unordered_map<string, string> read_fasta_file(const string &filename) {
  unordered_map<string, string> sequences;
  ifstream file(filename);

  if (file.is_open()) {
    string line;
    string current_key;

    while (getline(file, line)) {
      if (line[0] == '>') {
        current_key = line.substr(1);
      } else {
        sequences[current_key] += line;
      }
    }

    file.close();
  }

  return sequences;
}

string reverse_complement(const string &seq) {
  string rev_complement;
  rev_complement.reserve(seq.size());

  for (auto it = seq.rbegin(); it != seq.rend(); ++it) {
    switch (*it) {
      case 'A':
        rev_complement.push_back('T');
        break;
      case 'T':
        rev_complement.push_back('A');
        break;
      case 'C':
        rev_complement.push_back('G');
        break;
      case 'G':
        rev_complement.push_back('C');
        break;
      default:
        break;
    }
  }

  return rev_complement;
}


void save_predictions(const string &filename, const vector<string> &predictions) {
  ofstream file(filename);

  if (file.is_open()) {
    for (const auto &pred : predictions) {
      file << pred << endl;
    }
    file.close();
  }
}


int main() {
    string bound_file = "bound.fasta";
    string notbound_file = "notbound.fasta";
    string test_file = "test.fasta";

    auto bound_sequences_map = read_fasta_file(bound_file);
    auto notbound_sequences = read_fasta_file(notbound_file);
    auto test_sequences = read_fasta_file(test_file);


    vector<string> bound_sequences;
      for (const auto &seq_pair : bound_sequences_map) {
        bound_sequences.push_back(seq_pair.second);
      }


    int k = 15;
    vector<vector<double>> pwm = greedy_motif_search(bound_sequences, k);

    unordered_map<string, string> rev_complements;
    for (const auto &seq_pair : test_sequences) {
      rev_complements[seq_pair.first] = reverse_complement(seq_pair.second);
    }
    
    vector<pair<string, double>> scores;
    for (const auto &seq_pair : test_sequences) {
      double score = calculate_pwm_score(seq_pair.second, pwm);
      double rev_complement_score = calculate_pwm_score(rev_complements[seq_pair.first], pwm);
      scores.push_back({seq_pair.first, max(score, rev_complement_score)});
    }


    sort(scores.begin(), scores.end(), [](const auto &a, const auto &b) {
      return a.second > b.second;
    });

    vector<string> top_predictions;
    for (size_t i = 0; i < 2000 && i < scores.size(); ++i) {
      top_predictions.push_back(scores[i].first);
    }


    string output_file = "predictions.txt";
    save_predictions(output_file, top_predictions);


  return 0;
}
