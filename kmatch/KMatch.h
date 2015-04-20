#ifndef KMATCH_INCLUDED
#define KMATCH_INCLUDED 1
#include <vector>
#include <algorithm>
#include <iostream>
#include <fstream>
#define __STDC_LIMIT_MACROS
#include <stdint.h>
#include <cstdlib>

#define KMATCH_MAX_FREQ 4
#define KMATCH_NUC_A 0
#define KMATCH_NUC_C 1
#define KMATCH_NUC_G 2
#define KMATCH_NUC_T 3
#define KMATCH_NOKMER INT64_MIN

typedef struct kmer_positions_s {
  uint64_t kmer; //canonical kmer
  uint16_t seq_index;
  int32_t position; //sign indicates direction +FW -REV
  bool operator <(const kmer_positions_s& rhs) const{
    return (kmer < rhs.kmer);
  }

} kmer_positions_t;

typedef struct {
  int64_t t_position;
  int64_t q_position;
} kmer_matches_t; // if a kmer is High Frequency, then a (pos, -1) and a (-1, pos) are added. TODO: allow for HF in only one.

typedef struct {
  std::string name;
  uint64_t length;
} seq_attributes_t;

inline int64_t str_to_kmer(const char * _str); //returns a canonical kmer with sign indicating position, KMATCH_NOKMER if invalid input

class KMatch {
  public:
    KMatch(char * _target_filename, char * _query_filename, uint8_t _K);
    void load_positions();//loads both kmer arrays
    void merge_positions();
  private:
    uint8_t K;
    char * target_filename;
    char * query_filename;
    void kmer_array_from_fasta(char * filename, std::vector<kmer_positions_t> & kposv, std::vector<seq_attributes_t> & seqnames);
    
    std::vector<seq_attributes_t> target_seqs, query_seqs;
    std::vector<kmer_positions_t> target_positions, query_positions;
};

#endif //KMATCH_INCLUDED
