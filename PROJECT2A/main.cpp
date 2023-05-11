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
    string pwm_file = "project2a_PWM.txt";

    auto bound_sequences = read_fasta_file(bound_file);
    auto notbound_sequences = read_fasta_file(notbound_file);
    auto test_sequences = read_fasta_file(test_file);
    auto pwm = read_pwm_file(pwm_file);

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
