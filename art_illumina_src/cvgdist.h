#ifndef CVGDIST_H
#define CVGDIST_H

#include <map>
#include <string>
#include <vector>
using namespace std;

class cvgdist {
public:
  cvgdist();
  // read target coverage profile from file
  // FASTA-style format with '#' delimiting numbers, e.g.:
  //   >seq1
  //   1#2#5#12#10#10#[...]2#1#
  bool setdist(char* file_name);
  // check if sequence is present in profile
  bool has_seq(const string seq_id);
  // specify the sequence that is currently being sampled from
  void setseq(const string seq_id);
  // get random position within sequence using target coverage as weights
  long get_ran_pos_weighted();
  // get random position within sequence using target coverage as weights
  long get_ran_pos_weighted(const string &seq_id);
  ~cvgdist();

protected:
  // pick a random key from a map using values as cumulative probabilities (i.e. increase up to 1!)
  int get_ran_key_weighted(const map<int, double> &key2wgt);
  void debug_info();
  void debug_test();
  // stores for each sequence a map of copy number -> positions
  map< string, map< int, vector<long> > > seq2cvg2pos;
  // stores for each sequence a map of copy number -> weight
  map< string, map< int, double> > seq2cvg2wgt;
  // the sequence that is currently being sampled from
  string current_seq;
};

#endif /* CVGDIST_H */
