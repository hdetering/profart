#include <iostream>
#include <math.h>        /* floor */
#include <stdlib.h>      /* atoi */
#include "cvgdist.h"
#include "art.h"
#include "readSeqFile.h"

cvgdist::cvgdist() {}
cvgdist::~cvgdist() {}

bool cvgdist::setdist(char* file_name) {
  readSeqFile seq_reader(file_name);
  string id;
  string seq;

  while (seq_reader.next_seq(id, seq)) {
    long pos = 0;
    map<int, vector<long> > map_cvg2pos;
    string s_cvg;
    for(string::iterator it=seq.begin(); it!=seq.end(); ++it) {
      if (*it == '#') { // coverage values are separated by '#'
        int cvg = atoi(s_cvg.c_str());
        map<int, vector<long> >::iterator it_c2p = map_cvg2pos.find(cvg);
        if (it_c2p == map_cvg2pos.end()) { map_cvg2pos[cvg] = vector<long>(); }
        map_cvg2pos[cvg].push_back(pos);

        ++pos;
        s_cvg = "";
      }
      else {
        s_cvg += *it;
        continue;
      }
    }

    // calculate weights as fraction of sequence covered at X times X
    map<int, double> map_cvg2wgt;
    typedef map<int, vector<long> >::iterator it_type;
    double cumsum = 0.0;
    for (it_type it=map_cvg2pos.begin(); it!=map_cvg2pos.end(); ++it) {
      double w = it->first * (it->second.size() / double(pos));
      map_cvg2wgt[it->first] = w;
      cumsum += w;
    }
    // normalize and sum weights so they increase up to one
    double cp = 0.0;
    typedef map<int, double>::iterator it_mid;
    for (it_mid it=map_cvg2wgt.begin(); it!=map_cvg2wgt.end(); ++it) {
      cp += ( (it->second) / cumsum );
      it->second = cp;
    }

    this->seq2cvg2pos[id] = map_cvg2pos;
    this->seq2cvg2wgt[id] = map_cvg2wgt;
  }
debug_info();
//debug_test();
  return true;
}

bool cvgdist::has_seq(const string seq_id) {
  bool in_seq2cvg2pos = (this->seq2cvg2pos.find(seq_id) != this->seq2cvg2pos.end());
  bool in_seq2cvg2wgt = (this->seq2cvg2wgt.find(seq_id) != this->seq2cvg2wgt.end());
  return ( in_seq2cvg2pos && in_seq2cvg2wgt );
}

// set current sequence to process, make sure sequence exists
void cvgdist::setseq(const string seq_name) {
  map<string, map<int, double> >::iterator it = this->seq2cvg2wgt.find(seq_name);
  if (it != this->seq2cvg2wgt.end()) {
    this->current_seq = seq_name;
  }
  else {
    cerr << "[ERROR] sequence id '" << seq_name << "' not found in coverage distribution.";
    exit(1);
  }
}

long cvgdist::get_ran_pos_weighted() {
  return get_ran_pos_weighted(this->current_seq);
}

long cvgdist::get_ran_pos_weighted(const string &seq_name) {
  // pick coverage bucket
  int cvg = get_ran_key_weighted(this->seq2cvg2wgt[seq_name]);
  // pick position within coverage bucket
  int idx_ran = floor( r_prob() * this->seq2cvg2pos[seq_name][cvg].size() );

  return this->seq2cvg2pos[seq_name][cvg][idx_ran];
}

// pick a random key from a map using values as cumulative probabilities (i.e. uncrease up to 1!)
int cvgdist::get_ran_key_weighted(const map<int, double> &key2wgt) {
  double r = r_prob();
  map<int, double>::const_iterator it = key2wgt.begin();
  while (r > it->second && it != key2wgt.end()) { it++; }
  return it->first;
}

void cvgdist::debug_info() {
  typedef map< string, map< int, vector<long> > >::iterator it_m_sil;
  typedef map<int, vector<long> >::iterator it_m_il;
  cout << "Loaded coverage profile." << endl;
  cout << "  sequences: " << this->seq2cvg2pos.size() << endl;

  for (it_m_sil it1=this->seq2cvg2pos.begin(); it1!=this->seq2cvg2pos.end(); ++it1) {
    int i=0;
    double p=0.0;
    cout << "  " << it1->first << endl;
    for (it_m_il it2=it1->second.begin(); it2!=it1->second.end(); ++it2) {
      cout << "    " << it2->first << "\t" << it2->second.size() << "\t(" << this->seq2cvg2wgt[it1->first][it2->first]-p << ")" << endl;
      ++i;
      p = this->seq2cvg2wgt[it1->first][it2->first];
    }
  }
  cout << "" << endl;
}

void cvgdist::debug_test() {
  // sample some random positions and output ratios
  int rep = 100000; // replicates per sequence

  cout << "Picking random positions according to coverage weights..." << endl;
  typedef map< string, map< int, vector<long> > >::iterator it_m_sil;
  for (it_m_sil it=this->seq2cvg2pos.begin(); it!=this->seq2cvg2pos.end(); ++it) {
    cout << "  " << it->first << endl;
    map<int, int> map_pos2cvg; // remember hit positions
    for (int i=0; i<rep; ++i) {
      int pos = get_ran_pos_weighted(it->first);
      map<int, int>::iterator it_p2c = map_pos2cvg.find(pos);
      if (it_p2c == map_pos2cvg.end()) { map_pos2cvg[pos] = 0; }
      map_pos2cvg[pos]++;
    }
    for (map<int,int>::iterator it=map_pos2cvg.begin(); it!=map_pos2cvg.end(); ++it) {
      cout << "    " << it->first << ":\t" << it->second/double(rep) << endl;
    }
  }
}
