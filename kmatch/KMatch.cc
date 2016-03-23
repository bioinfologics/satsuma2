#include "KMatch.h"
#include <sys/time.h>
#include <thread>
//#include <functional>

void timed_log(std::string s){
  struct timeval tp;
  gettimeofday(&tp,NULL);
  long int ms = tp.tv_sec * 1000 + tp.tv_usec / 1000;
  std::cout<<"TIME_LOG: "<<ms<<" - "<<s<<std::endl;
};

KMatch::KMatch(char * _target_filename, char * _query_filename, uint8_t _K, int _max_freq){
  target_filename=_target_filename;
  query_filename=_query_filename;
  K=_K;
  max_freq=_max_freq;
};

void KMatch::kmer_array_from_fasta(char * filename, std::vector<kmer_position_t> & kposv, std::vector<seq_attributes_t> & seqnames){
  std::vector<kmer_position_t> karray;
  //read fasta and push_back(kmer,pos) (use pos as chr*CHR_CONST+offset)
  //open file
  std::string line,seq;
  std::ifstream fasta(filename);
  seq_attributes_t seq_attr;
  timed_log(" Loading array ");
  std::pair<uint64_t,bool> ckmer;
  int64_t chr_offset=0;
  uint32_t seq_index=0;
  uint64_t kmer_index=0;
  const uint64_t KMER_MASK=( ((uint64_t)1)<<(K*2) )-1;
  const uint64_t KMER_FIRSTOFFSET=(K-1)*2;
  std::cout<<"Loading fasta "<<filename<<" into kmer array"<<std::endl;
  seq="";
  while ( 1 ){
    getline (fasta, line);
    if ( (line.size()>0 && line[0]=='>') || fasta.eof()){
      if (seq.size()>0) {//There is a previous sequence, process it!
        //std::cout<<"processing sequence"<<std::endl;
        seq_attr.length=seq.size();
        kposv.resize(kposv.size()+seq.size()+1-K);//XXX: this could be optimised to at least grow N positions if growth needed, so it doesn't grow in every small sequence
        seqnames.push_back(seq_attr);
        const char * s=seq.c_str();
        //TODO: further speedup? reserve all space first, choose a kmer value as threshold and insert larger kmers fromt the top, smaller or equal from the bottom, sort the two parts (frontier will be the next-insertion point) independently and join them in the filtering step.
        int64_t last_unknown=-1;
        uint64_t fkmer=0,rkmer=0;
        for (uint64_t p=0;p<seq.size();p++){
          //fkmer: grows from the right (LSB) 
          //rkmer: grows from the left (MSB)
          switch (s[p]) {
            case 'A':
            case 'a':
              fkmer=( (fkmer<<2) + 0 ) & KMER_MASK;
              rkmer=(rkmer>>2) + (((uint64_t) 3)<<KMER_FIRSTOFFSET) ;
              break;
            case 'C':
            case 'c':
              fkmer=( (fkmer<<2) + 1 ) & KMER_MASK;
              rkmer=(rkmer>>2) + (((uint64_t) 2)<<KMER_FIRSTOFFSET) ;
              break;
            case 'G':
            case 'g':
              fkmer=( (fkmer<<2) + 2 ) & KMER_MASK;
              rkmer=(rkmer>>2) + (((uint64_t) 1)<<KMER_FIRSTOFFSET) ;
              break;
            case 'T':
            case 't':
              fkmer=( (fkmer<<2) + 3 ) & KMER_MASK;
              rkmer=(rkmer>>2) + (((uint64_t) 0)<<KMER_FIRSTOFFSET) ;
              break;
            default:
              fkmer=( (fkmer<<2) + 0 ) & KMER_MASK;
              rkmer=(rkmer>>2) + (((uint64_t) 3)<<KMER_FIRSTOFFSET);
              last_unknown=p;
              break;
          }
          //std::cout<<"c="<<s[p]<<" f="<<fkmer<<" r="<<rkmer<<std::endl;
          //TODO: last unknown passed by?
          if (last_unknown+K<=p){
            /*char cstring[32],fkstring[32],rkstring[32];
            cstring[K]=0;
            fkstring[K]=0;
            rkstring[K]=0;
            for (int i=0;i<K;i++){
              cstring[i]=s[p-K+i+1];
              fkstring[i]='0'+( fkmer >> (2*(K-i-1)) )%4;
              rkstring[i]='0'+( rkmer >> (2*(K-i-1)) )%4;
            };
            std::cout<<"string="<<cstring<<" fkmer="<<fkmer<<" ("<<fkstring<<") rkmer="<<rkmer<<" ("<<rkstring<<")"<<std::endl;*/
            //result is min(kmer/rkmer), and set position / reverse
            if (fkmer<=rkmer){
              kposv[kmer_index].kmer=fkmer;
              kposv[kmer_index].position=p-K+2+seq_index*KMATCH_POSITION_CHR_CNST;//1-based position
            } else {
              kposv[kmer_index].kmer=rkmer;
              kposv[kmer_index].position=-(p-K+2+seq_index*KMATCH_POSITION_CHR_CNST);//1-based position
            }
          kmer_index++;
          }/*else if (p>=K){ XXX:review this and make it ok for the start!
            //result is NOKMER
            //kposv[kmer_index].kmer=KMATCH_NOKMER;
            //kposv[kmer_index].position=kmer_index+1;//1-based position

          }*/

        }
      }
      if (fasta.eof()){
        break;
      } else {//init new seq_attributes;
        //std::cout<<"Loading sequence '"<<line<<"'"<<std::endl;
        seq_attr.name=line;
        seq="";
        seq_index++;
      }
    }
    else {
      seq+=line;
    }
  }
  fasta.close();
  std::cout<<"Kmer array with "<<kmer_index<<" elements created"<<std::endl;
  timed_log(" Sorting array ");
  std::sort(kposv.begin(),kposv.end());
  std::cout<<"Kmer array sorted"<<std::endl;
  //TODO: allow only K elements, replace with KMATCH_NOKMER if an element is too high frequency?
  uint64_t ri,wi=0,f;
  uint64_t kposv_size=kposv.size();
  for (ri=0;ri<kposv_size;ri++){
    for (f=0;ri+f<kposv_size && kposv[ri+f].kmer==kposv[ri].kmer;f++);
    if (f>max_freq) {
      ri+=f-1;//jump to the last one
      //std::cout<<"filtering equal elements ("<<kposv[ri].kmer<<") at "<<ri<<"-"<<ri+f-1<<std::endl;
    } else {
      kposv[wi++]=kposv[ri];
    }
  }
  kposv.resize(wi-1);
  std::cout<<"Kmer array filtered to "<<wi<<" elements"<<std::endl;
}

void KMatch::load_target_positions(){
  kmer_array_from_fasta(target_filename,target_positions,target_seqs);
}
void KMatch::load_query_positions(){
  kmer_array_from_fasta(query_filename,query_positions,query_seqs);
}

void KMatch::merge_positions(){
  //TODO: parallel, each thread can take 1/t of the first array and scan the second array until finding the suitable start.
  std::cout<<"Starting to create matching positions"<<std::endl;
  uint64_t ti=0;
  kmer_match_t m;
  uint64_t tsize=target_positions.size();
  uint64_t qsize=query_positions.size();
  for (uint64_t qi=0;qi<qsize;qi++){
    //advance target
    while (ti<tsize && query_positions[qi].kmer>target_positions[ti].kmer) ti++;
    //match on target?
    if (query_positions[qi].kmer==target_positions[ti].kmer){
      //check for multi-match
      for (uint64_t j=0; ti+j<tsize && query_positions[qi].kmer==target_positions[ti+j].kmer; j++) {
        bool rev=false;
        //insert each result XXX insert the real position + 1 to allow for sign
        if (query_positions[qi].position>0){
          m.q_position=query_positions[qi].position;
        }else{
          rev=true;
          m.q_position=-query_positions[qi].position;
        }
        if (target_positions[ti+j].position>0){
          m.t_position=target_positions[ti+j].position;
        } else {
          rev=!rev;
          m.t_position=-target_positions[ti+j].position;
        }
        m.reverse=rev;
        kmatches.push_back(m);//XXX: should this be pre-allocated?
      }
    }
  }
  std::cout<<kmatches.size()<<" matching positions"<<std::endl;

}

void KMatch::clear_positions(){
  target_positions.clear();
  query_positions.clear();
}

void KMatch::dump_matching_blocks(char * out_filename, int min_length, int max_jump){
  //XXX: allow for multi matches!!
  //watch out: matching positions are 1-based to allow for sign always
  uint64_t kmsize=kmatches.size();
  std::cout<<"Sorting the matches"<<std::endl;
  std::sort(kmatches.begin(),kmatches.end());//XXX: this needs to be sorted by absolute value!
  std::cout<<"Dumping matches of "<<min_length-K+1 << " kmers with jumps of up to "<<max_jump<<"kmers to "<<out_filename<<std::endl;
  std::ofstream out_file(out_filename);
  int64_t match_start=0;
  int64_t q_delta=0, t_delta=1;//Invalid values, just to make sure.
  uint64_t dumped=0;
  std::list<multikmer_match_t> active_matches;
  //MULTI-MATCH: keep a list of "started matches" with q_start, t_start,  last-match
  for (uint64_t i=1;i<=kmsize ;i++){//do not check the last element, check it outside!
    //MULTI-MATCH version:
    bool used=0;
    //1) for each match in "current matches":
    //std::cout<<"--\nEvaluating kmatch on "<< kmatches[i].q_position <<"->"<< kmatches[i].t_position <<" (r="<<kmatches[i].reverse<<")"<<std::endl;
    for (auto am=active_matches.begin();am!=active_matches.end();){
      //std::cout<<" Existing match on "<< am->q_start <<"->"<< am->t_start <<" (r="<<am->reverse<<")"<<std::endl;
      q_delta=kmatches[i].q_position-am->q_start;
      t_delta=(am->reverse ? am->t_start - kmatches[i].t_position : kmatches[i].t_position-am->t_start);
      //std::cout<<" qd="<<q_delta<<" td="<<t_delta<<std::endl;
      
      //  b) if current kmer-pair extends, update last-match
      if (q_delta == t_delta && q_delta-am->length <=max_jump && //check coordinates
          am->reverse == kmatches[i].reverse ) { //check orientation
          am->length=q_delta;
          used=1;
          //std::cout<<"Match updated!!!"<<std::endl;
      }
        
      //  a) if current kmer-pair[0] out of jump range, check (dump) and remove from list
      if (q_delta - am->length >max_jump || i==kmsize) {//match has finished!
        //std::cout<<" Match is past its jump "<<std::endl;
        if (am->length+K >= min_length ){ //check min_length
          //std::cout<<" Dumping match (length="<<am->length+K<< ")"<<std::endl;
          //dump!
          t_result r;
          r.query_id=am->q_start/KMATCH_POSITION_CHR_CNST-1;//0-based
          r.target_id=am->t_start/KMATCH_POSITION_CHR_CNST-1;//0-based
          r.query_size=query_seqs[r.query_id].length;
          //XXX:review position and K displacement when reverse, etc
          r.qstart=am->q_start%KMATCH_POSITION_CHR_CNST-1;//0-based
          if (am->reverse) {
            r.tstart=(am->t_start-am->length)%KMATCH_POSITION_CHR_CNST-1;//0-based
          } else {
            r.tstart=am->t_start%KMATCH_POSITION_CHR_CNST-1;//0-based
          }
          r.len=am->length+K;
          r.reverse=am->reverse;
          r.prob=1;
          r.ident=1;
          out_file.write((char *) &r,sizeof(r));
          dumped++;
        }
        //delete from matches
        am=active_matches.erase(am);
	//if (am!=active_matches.begin()) am--;
      } else am++;
    }
    //2) If current kmer didnt extend any match,start a new match with it
    if (!used){
      multikmer_match_t mkm;
      mkm.q_start=kmatches[i].q_position;
      mkm.t_start=kmatches[i].t_position;
      mkm.length=0;
      mkm.reverse=kmatches[i].reverse;
      active_matches.push_back(mkm);
      //std::cout<<"NEW match on "<< mkm.q_start <<"->"<<mkm.t_start<<std::endl;
    }
    




    //TODO: to allow multi-matches just change i-1 for a back-search of the previous link in this match (i.e, go back till position[j] <position-max_jmp and if any point matches, use it to move forward.
    //TODO: to allow for multi-matches the start register needs to change to a vector with some more variables.
    //if match breaks in this element
    /*----------------
    if (i<kmsize){
      q_delta=kmatches[i].q_position-kmatches[i-1].q_position;
      t_delta=kmatches[i].t_position-kmatches[i-1].t_position;
    }
    if ( i==kmsize //end condition!
         || kmatches[i].reverse != kmatches[i-1].reverse //change in direction
         || q_delta-1>max_jump  //long jump
         || (kmatches[i].reverse==false && q_delta != t_delta)//direct and not same difference
         || (kmatches[i].reverse==true && q_delta != -t_delta) ) { //reverse and not same difference 

      //std::cout<<"evaluating match in ["<<match_start<<"-"<<i-1<<"] "<<q_delta<<" "<<t_delta<<" "<<kmatches[i].reverse<<" "<<kmatches[i-1].reverse<<"||"<<(kmatches[i].reverse != kmatches[i-1].reverse)<<" "<<(q_delta-1>max_jump)<<" "<< (kmatches[i].reverse==false && q_delta != t_delta)<<" "<<(kmatches[i].reverse==true && q_delta != -t_delta)<<std::endl;
      //length>min_length?
      if (i-match_start>=min_length-K){
        //TODO:dump
        t_result r;
        r.query_id=kmatches[match_start].q_position/KMATCH_POSITION_CHR_CNST-1;//0-based
        r.target_id=kmatches[match_start].t_position/KMATCH_POSITION_CHR_CNST-1;//0-based
        r.query_size=query_seqs[r.query_id].length;
        //XXX:review position and K displacement when reverse, etc
        r.qstart=kmatches[match_start].q_position%KMATCH_POSITION_CHR_CNST-1;//0-based
        if (kmatches[match_start].reverse) {
          r.tstart=kmatches[i-1].t_position%KMATCH_POSITION_CHR_CNST-1;//0-based
        } else {
          r.tstart=kmatches[match_start].t_position%KMATCH_POSITION_CHR_CNST-1;//0-based
        }
        r.len=kmatches[i-1].q_position - kmatches[match_start].q_position+K;
        r.reverse=kmatches[match_start].reverse;
        r.prob=1;
        r.ident=1;
        out_file.write((char *) &r,sizeof(r));
        dumped++;
        //std::cout<<"Match dumped! "<<r.query_id<<"["<<r.qstart<<":"<<r.qstart+r.len<<"] -> "<<r.target_id<<"["<<r.tstart<<":"<<r.tstart+r.len<<"]"<<std::endl;
      }
      //move match start
      match_start=i;
    }
  ----------*/
  }
  out_file.close();
  std::cout<<dumped<<" matches dumped"<<std::endl;
}

int main(int argc, char ** argv){
  if (argc!=8) {
    std::cout<<"Usage: "<<argv[0]<<" query.fa target.fa K output.fa min_length jump max_freq"<<std::endl;
    return -1;
  }
  
  if (atoi(argv[3])%2==0) {
    std::cout<<"KMatch only accepts odd K values, please try again"<<std::endl;
    return 1;
  }
  //timed_log(" START ");
  KMatch kmatch(argv[2],argv[1],atoi(argv[3]),atoi(argv[7]));
  //timed_log(" load_positions() ");
  std::thread q(&KMatch::load_query_positions,std::ref(kmatch));
  std::thread t(&KMatch::load_target_positions,std::ref(kmatch));
  q.join();
  t.join();
  //timed_log(" merge_positions() ");
  kmatch.merge_positions();
  //timed_log(" clear_positions() ");
  kmatch.clear_positions();
  //timed_log(" dump_matching_blocks() ");
  kmatch.dump_matching_blocks(argv[4],atoi(argv[5]),atoi(argv[6]));
  //timed_log(" END ");
}
