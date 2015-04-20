#include "kmatch/KMatch.h"

inline int64_t str_to_kmer(const char * _str,uint8_t _K){
  int64_t key=0,rkey=0;
  for (uint8_t i=0;i<_K;i++){
    key<<=2;
    switch(_str[i]){
      case 'A':
      case 'a':
        key+=KMATCH_NUC_A;
        rkey+=((int64_t) KMATCH_NUC_T)<<(2*i);
        break;
      case 'C':
      case 'c':
        key+=KMATCH_NUC_C;
        rkey+=((int64_t) KMATCH_NUC_G)<<(2*i);
        break;
      case 'G':
      case 'g':
        key+=KMATCH_NUC_G;
        rkey+=((int64_t) KMATCH_NUC_C)<<(2*i);
        break;
      case 'T':
      case 't':
        key+=KMATCH_NUC_T;
        rkey+=((int64_t) KMATCH_NUC_A)<<(2*i);
        break;
      default:
        //XXX: not fancy at all!!!
        return KMATCH_NOKMER;
        break;
    }
  }
  return (key < rkey ? key : -rkey);

};

KMatch::KMatch(char * _target_filename, char * _query_filename, uint8_t _K){
  target_filename=_target_filename;
  query_filename=_query_filename;
  K=_K;
};

void KMatch::kmer_array_from_fasta(char * filename, std::vector<kmer_positions_t> & kposv, std::vector<seq_attributes_t> & seqnames){
  std::vector<kmer_positions_t> karray;
  //read fasta and push_back(kmer,pos) (use pos as chr*CHR_CONST+offset)
  //open file
  std::string line,seq;
  std::ifstream fasta(filename);
  seq_attributes_t seq_attr;
  int64_t kmer;
  int64_t chr_offset=0;
  uint32_t seq_index;
  uint64_t kmer_index=0;
  //while (!EOF)
  std::cout<<"Loading fasta "<<filename<<" into kmer array"<<std::endl;
  while ( getline (fasta, line)){
    if ( (line.size()>0 && line[0]=='>') || fasta.eof()){
        if (seq.size()>0) {seq_attr.length=seq.size();
        kposv.resize(kposv.size()+seq.size()+1-K);
        seqnames.push_back(seq_attr);
        const char * s=seq.c_str();
        for (uint64_t p=0;p<seq.size()+1-K;p++){
          //TODO: this could well be parallel
          kposv[kmer_index].seq_index=seq_index;
          kposv[kmer_index].kmer=str_to_kmer(s+p,K);
          kposv[kmer_index].position=p;
          kmer_index++;
        }
      }
      if (line.size()>0 && line[0]=='>'){
        //init new seq_attributes;
        std::cout<<"Loading sequence '"<<line<<"'"<<std::endl;
        seq_attr.name=line;
        seq="";
        seq_index++;
      }
    }
    else {
      seq+=line;
    }
  }
  std::cout<<"Kmer array with "<<kmer_index<<" elements created"<<std::endl;
  std::sort(kposv.begin(),kposv.end());
  std::cout<<"Kmer array sorted"<<std::endl;
}

void KMatch::load_positions(){
  kmer_array_from_fasta(target_filename,target_positions,target_seqs);
  kmer_array_from_fasta(query_filename,query_positions,query_seqs);
}

int main(int argc, char ** argv){
  KMatch kmatch(argv[1],argv[2],atoi(argv[3]));
  kmatch.load_positions();
}
