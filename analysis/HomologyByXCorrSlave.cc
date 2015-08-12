#ifndef FORCE_DEBUG
#define NDEBUG
#endif

#include "analysis/DNAVector.h"
#include <string>
#include <unistd.h>
#include "base/CommandLineParser.h"
#include "analysis/CrossCorr.h"
#include "analysis/SequenceMatch.h"
#include "analysis/SeqChunk.h"
#include "analysis/AlignProbability.h"
#include "analysis/MatchDynProg.h"
#include "analysis/WorkQueue.h"
#include "analysis/ProbTable.h"
#include "util/mutil.h"
#include <queue>
#include <thread>
#include <mutex>
#include <sys/types.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <netdb.h>

//global variables, cheating but should work

/*FILE input is done once */
vecDNAVector targetRaw, queryRaw; //holds the sequences
svec<string> targetNames, queryNames; //holds the sequences names
ChunkManager * cmQuery; //XXX: forced to use this as no default init
ChunkManager * cmTarget; //XXX: forced to use this as no default init
vecDNAVector target, query;
svec<SeqChunk> targetInfo, queryInfo;
double targetTotal;
unsigned int debug_max_targets;

std::queue<t_pair> targets;
std::mutex targets_mutex;
std::mutex creation_mutex;
std::vector<t_result> results;
std::mutex results_mutex;
std::mutex chunking_mutex;
int minLen;
int targetChunk;
int queryChunk;
double minProb;

ProbTable probt;
double topCutoff;
double topCutoffFast;
bool prob_table;
bool processing_finished=false;
string target_filename,query_filename;

//====== XXX: REALLY????? =======================================
int RCQuery(int offset, int start, int len, int chunkSize, int qLen) 
{
  int l = start + qLen - offset - chunkSize;
  return l;
}

class HomologyByXCorr
{
  public:
    HomologyByXCorr() {
      m_minLen = minLen;
      m_pMulti = NULL;
      m_offset = 0;
      m_qLen = 0;
      m_chunk = 0;
      m_tOffset = 0;
      m_qID = -1;
      m_tID = -1;

      m_targetSize = 0;
      m_minProb = 0.99;
    }

    void SetMinProb(double p) {
      m_minProb = p;
    }

    void Align(int _target_id, const CCSignal & _target_signal, int _query_id, const CCSignal & _query_signal, const CCSignal & _query_rcsignal, bool fast);

    void SetTargetSize(double i) {
      m_targetSize = i;
    }

    void SetMinimumAlignLen(int l) {
      m_minLen = l;
    }

    void disconnect_from_master();
    int connect_to_master();
    int get_targets_from_master();
    void send_solutions_to_master();
    std::vector<CCSignal> create_signals(const vecDNAVector & _sequences, unsigned long int _from, unsigned long int _count, bool reverse);
    void align_target(t_pair p);
    void FilterMatches(int _target_id, int _query_id, const DNAVector & _query_seq, const vecSeqMatch & _matches, bool _reverse);
    void work();


  private:
    int sid;
    int threads;
    int m_minLen;
    svec<int> m_bestID;
    svec<double> m_ident;
    MultiMatches * m_pMulti;

    int m_tID;
    int m_qID;
    int m_tOffset;
    int m_offset;
    int m_qLen;
    int m_chunk;
    double m_targetSize;
    double m_minProb;

    int sockfd;
    string master_hostname;
    int portno;
    unsigned int slave_id;
    int num_targets;



};

void HomologyByXCorr::disconnect_from_master(){
  //TODO: should we clear the target-pair queue here?
  close(sockfd);
  sockfd=0;
}

int HomologyByXCorr::connect_to_master(){
  //connects to the master, returns sockfd or -1 for failure
  //TODO: sends its own id to master?
  int sockfd;
  //int  portno, n;
  struct sockaddr_in serv_addr;
  struct hostent *server;

  sockfd = socket(AF_INET, SOCK_STREAM, 0);
  if (sockfd <= 0){
    cout << "ERROR opening socket" << endl;
    return -1; //XXX returning error code
  }
  server = gethostbyname(master_hostname.c_str());
  if (server == NULL) {
    cout << "ERROR resolving master hostname" << endl;
    return -1;
  }
  bzero((char *) &serv_addr, sizeof(serv_addr));
  serv_addr.sin_family = AF_INET;
  bcopy((char *)server->h_addr, (char *)&serv_addr.sin_addr.s_addr, server->h_length);
  serv_addr.sin_port = htons(portno);
  if (connect(sockfd,(struct sockaddr *) &serv_addr, sizeof(serv_addr)) < 0) {
    cout << "ERROR connecting to master" << endl;
    close(sockfd);
    return -1;
  }
  //TODO: write own ID to server and wait for an OK

  return sockfd;

}
void HomologyByXCorr::FilterMatches(int _target_id, int _query_id, const DNAVector & _query_seq, const vecSeqMatch & _matches, bool _reverse){
  int len, tStart, qStart;
  std::vector<t_result> local_results;
  for (int j=0; j<_matches.isize(); j++) {
    if (_matches[j].GetLength() < m_minLen)
      continue;

    len = _matches[j].GetLength();
    tStart = targetInfo[_target_id].GetStart() + _matches[j].GetStartTarget();
    if (!_reverse) {
      qStart = queryInfo[_query_id].GetStart() +  _matches[j].GetStartQuery();
    }else{
      qStart = RCQuery(queryInfo[_query_id].GetStart(), _matches[j].GetStartQuery(), len, queryChunk, cmQuery->GetSize(queryInfo[_query_id].GetID()));
    }
    double prob;
    double ident;
    if (prob_table){
      prob=probt.GetMatchProbability( ident, target[_target_id], _query_seq, _matches[j].GetStartTarget(), _matches[j].GetStartQuery(), _matches[j].GetLength());
      if (prob < m_minProb) continue;
    } else {
      prob = GetMatchProbability(target[_target_id],
        _query_seq,
        _matches[j].GetStartTarget(),
        _matches[j].GetStartQuery(),
        _matches[j].GetLength(),
        targetTotal);
      if (prob < m_minProb) continue;
      ident = PrintMatch(_query_seq, target[_target_id], _matches[j], true);//XXX: does this really print?
    }



    //cout << "Match # " << j << " probability " << 100.* prob << " %" << endl;
    //cout << "Start target: " << tStart << " - " << tStart + len << endl;
    t_result r;
    r.query_id=queryInfo[_query_id].GetID();
    r.target_id=targetInfo[_target_id].GetID();
    r.query_size=cmQuery->GetSize(r.query_id);
    r.qstart=qStart;
    r.tstart=tStart;
    r.len=len;
    r.reverse=_reverse;
    r.prob=prob;
    r.ident=ident;
    local_results.push_back(r);
  }
  //cout<<std::this_thread::get_id()<<" inserted matches after filtering: "<<local_results.size()<<endl;
  results_mutex.lock();
  results.insert( results.end(), local_results.begin(), local_results.end() );
  results_mutex.unlock();

}

void HomologyByXCorr::Align(
    int _target_id,
    const CCSignal & _target_signal, 
    int _query_id,
    const CCSignal & _query_signal, 
    const CCSignal & _query_rcsignal, 
    bool _fast)
{
  CrossCorrelation xc;
  vecSeqMatch matches,revmatches;
  SeqAnalyzer sa;
  int i, j;
  svec<float> result,revresult; 
  /*if (_fast && prob_table) {
    cout<< "WARNING: fast cutoff required, but probability table is being used!"<<endl;
  }*/
  sa.SetTopCutoff((_fast ? topCutoffFast : topCutoff));

  xc.CrossCorrelate(result, _target_signal, _query_signal);
  matches.clear();
  sa.MatchUp(matches, query[_query_id], target[_target_id], result);
  //cout<<std::this_thread::get_id()<<"FW matches before filtering: "<<matches.size()<<endl;
  FilterMatches(_target_id, _query_id, query[_query_id], matches, false);

  xc.CrossCorrelate(revresult, _target_signal, _query_rcsignal);
  revmatches.clear();
  DNAVector rcquery=query[_query_id];
  rcquery.ReverseComplement();
  sa.MatchUp(revmatches, rcquery, target[_target_id], revresult);
  //cout<<std::this_thread::get_id()<<"RC matches before filtering: "<<revmatches.size()<<endl;
  FilterMatches(_target_id, _query_id, rcquery, revmatches, true);

}
std::vector<CCSignal> HomologyByXCorr::create_signals( const vecDNAVector & _sequences, unsigned long int _from, unsigned long int _count, bool reverse){
  std::vector<CCSignal> signals;
  signals.resize(_count);
  for (unsigned long int i=0;i<_count;i++){
    if (!reverse){
      signals[i].SetSequence(_sequences[_from+i], targetChunk *2);//XXX should chunkSize be a parameter???
    }else{
      DNAVector revseq=_sequences[_from+i];
      revseq.ReverseComplement();
      signals[i].SetSequence(revseq, targetChunk *2);//XXX should chunkSize be a parameter???
    }
  }
  return signals;

}

void HomologyByXCorr::align_target(t_pair p){

  //cout<<"Align block running (block: - t:"<<p.targetFrom<<"-"<<p.targetTo<<" q:"<<p.queryFrom<<"-"<<p.queryTo<<")"<<endl;
  //XXX: is this storage correct? create all signals for the target and query spaces
  std::vector<CCSignal> targetSignals;
  std::vector<CCSignal> querySignals;
  std::vector<CCSignal> querySignalsReverse;
  unsigned long int target_count, query_count;

  //Create all signals
  //XXX: is To included or not?
  target_count=p.targetTo+1-p.targetFrom;
  query_count=p.queryTo+1-p.queryFrom;
  //XXX what is the coordinates, how are them set?
  targetSignals = create_signals(target, p.targetFrom, target_count, false);
  querySignals = create_signals(query, p.queryFrom, query_count, false);
  querySignalsReverse = create_signals(query, p.queryFrom, query_count, true);

  //cout<<"Block "<<pair_id<<" signals created..."<<endl;

  for (unsigned long int qi=0; qi<query_count; qi++){
    //for each query
    for (unsigned long int ti=0; ti<target_count; ti++){
      //for each target
      Align(
          p.targetFrom+ti, targetSignals[ti],
          p.queryFrom+qi, querySignals[qi], querySignalsReverse[qi],
          p.fast);
    }
  }
}

void HomologyByXCorr::work(){
  //cout <<" I finished creating the chunks!!!!!"<<endl;
  //Work loop, gets targets from master and reports back
  while(!processing_finished){
    //locks targets
    targets_mutex.lock();
    //pops a target
    if (!targets.empty()){
      t_pair p=targets.front();
      targets.pop();
      targets_mutex.unlock();
      align_target(p);
    } else {
      targets_mutex.unlock();
      sleep(1);
    }
  }
}

void launch_worker(){
  HomologyByXCorr hbxc;
  cout<<"worker created, now to work!!!"<<endl;
  hbxc.work();
}


//================================================================
//================================================================
//================================================================
int main( int argc, char** argv ){
  //======= Parse arguments =======
  //TODO: how much of this can come form the master with the instructions?
  commandArg<string> mStringCmmd("-master","name of the submit host");
  commandArg<int> portIntCmmd("-port","port of the master to connect", MYPORT);
  commandArg<int> sidIntCmmd("-sid","slave_id");
  commandArg<int> threadsIntCmmd("-p","number of working processes (threads)",1);
  commandArg<string> aStringCmmd("-q","query fasta sequence");
  commandArg<string> bStringCmmd("-t","target fasta sequence");
  commandArg<int> lIntCmmd("-l","minimum alignment length", 0);
  commandArg<int> qChunkCmmd("-q_chunk","query chunk size", 4096);
  commandArg<int> tChunkCmmd("-t_chunk","target chunk size", 4096);
  commandArg<double> probCmmd("-min_prob","minimum probability to keep match", 0.9999);
  commandArg<double> cutoffCmmd("-cutoff","signal selection cutoff", 1.8);
  commandArg<double> cutoffFastCmmd("-cutoff_fast","signal selection cutoff (fast)", 2.9);
  commandArg<bool> probtableCmmd("-prob_table","lookup probability in table (faster, less accurate)", false);
  commandArg<int> debugTargetsCmmd("-debug_targets","number of targets to receive before exiting (for debug purposes mostly, don't use)",0);


  commandLineParser P(argc,argv);
  P.SetDescription("Compares two sequences via cross-correlation");
  P.registerArg(mStringCmmd);
  P.registerArg(portIntCmmd);
  P.registerArg(threadsIntCmmd);
  P.registerArg(sidIntCmmd);
  P.registerArg(aStringCmmd);
  P.registerArg(bStringCmmd);
  P.registerArg(lIntCmmd);

  P.registerArg(qChunkCmmd);
  P.registerArg(tChunkCmmd);
  P.registerArg(probCmmd);
  P.registerArg(cutoffCmmd);
  P.registerArg(cutoffFastCmmd);
  P.registerArg(probtableCmmd);
  P.registerArg(debugTargetsCmmd);

  P.parse();

  string master_hostname = P.GetStringValueFor(mStringCmmd);
  int portno = P.GetIntValueFor(portIntCmmd);
  int slave_id = P.GetIntValueFor(sidIntCmmd);
  int threads = P.GetIntValueFor(threadsIntCmmd);
  query_filename = P.GetStringValueFor(aStringCmmd);
  target_filename = P.GetStringValueFor(bStringCmmd);
  minLen = P.GetIntValueFor(lIntCmmd);
  targetChunk = P.GetIntValueFor(tChunkCmmd);
  queryChunk = P.GetIntValueFor(qChunkCmmd);
  minProb = P.GetDoubleValueFor(probCmmd);

  topCutoff = P.GetDoubleValueFor(cutoffCmmd);
  topCutoffFast = P.GetDoubleValueFor(cutoffFastCmmd);
  prob_table = P.GetBoolValueFor(probtableCmmd);
  debug_max_targets = P.GetIntValueFor(debugTargetsCmmd);

  //======= Pre-Loading fasta files  =======
  time_t time_start = time(NULL);
  cout << "Loading query sequence:  " << query_filename << endl;
  queryRaw.Read(query_filename,queryNames);
  cout << " - Creating query chunks..."<<endl;
  cmQuery = new ChunkManager(queryChunk,0);
  cmQuery->ChunkIt(query, queryInfo, queryRaw, queryNames, 0, 0);
  queryRaw.clear();
  cout<<"DONE"<<endl;


  cout << "Loading target sequence:  " << target_filename << endl;
  targetRaw.Read(target_filename,targetNames);
  cout << " - Creating target chunks..."<<endl;
  cmTarget = new ChunkManager(targetChunk, targetChunk / 4);
  cmTarget->ChunkIt(target, targetInfo, targetRaw, targetNames, 0, 0);
  targetRaw.clear();

  //TODO: set total target size and check if needed to create the multi!!!!!!
  targetTotal = 0;
  for (int i=0; i<cmTarget->GetCount(); i++) {
    targetTotal += (double)cmTarget->GetSize(i);
  }
  cout<<"DONE"<<endl;
  cout<<"TIME SPENT ON LOADING: "<<(time(NULL) - time_start)<<endl;
  time_start = time(NULL);
  
  if (prob_table) {
    cout << " Initializing Probability table... for size "<< targetTotal<<" and cutoff "<< minProb <<endl;
    probt=ProbTable(targetTotal,minProb);
    //cout << " Initializing Probability table... for size "<< targetChunk<<" and cutoff "<< topCutoff <<endl;
    //probt=ProbTable(targetChunk,topCutoff);
    cout << "Done!!" << endl;
    cout<<"TIME SPENT ON CREATING PROBTABLE: "<<(time(NULL) - time_start)<<endl;
    time_start = time(NULL);
  }

  //======= Main loop
  cout<< "== launching workers =="<<endl;

  std::thread workers[threads];
  /*create_chunks(query_filename,target_filename);*/
  for (int wt=0;wt<threads;wt++) {
    workers[wt]= std::thread(launch_worker);
  }

  cout<< "== Entering communication loop =="<<endl;
  //server name resolution and such
  struct hostent *server;
  struct sockaddr_in serv_addr;
  cout << "comm loop for "<<master_hostname.c_str()<<" "<<portno<<endl;
  server = gethostbyname(master_hostname.c_str());
  if (server == NULL) {
    cout << "ERROR resolving master hostname" << endl;
    return -1;
  }
  bzero((char *) &serv_addr, sizeof(serv_addr));
  serv_addr.sin_family = AF_INET;
  bcopy((char *)server->h_addr, (char *)&serv_addr.sin_addr.s_addr, server->h_length);
  serv_addr.sin_port = htons(portno);
  unsigned int total_targets;
  while(!processing_finished){
    //acquire lock on target_blocks
    targets_mutex.lock();
    //check targets<threads*4?
    if (targets.size()<2*threads){
      targets_mutex.unlock();
      //PHASE 1: OPEN CONNECTION
      int sockfd=0;
      while (sockfd<=0) {
        sockfd = socket(AF_INET, SOCK_STREAM, 0);
        if (sockfd <= 0){
          cout << "ERROR opening socket" << endl;
          sleep(1);
          continue;
        }

        if (connect(sockfd,(struct sockaddr *) &serv_addr, sizeof(serv_addr)) < 0) {
          cout << "ERROR connecting to master: " <<strerror(errno)<< endl;
          close(sockfd);
          sockfd=0;
          sleep(1);
          continue;
        }
        if (write(sockfd,&slave_id,sizeof(slave_id)) != sizeof(slave_id)) {
          cout << "ERROR sending id to master" << endl;
          close(sockfd);
          sockfd=0;
          sleep(1);
          continue;
        }

      }
      //PHASE 2: SEND RESULTS
      //lock results
      results_mutex.lock();
      //create a single big-mem-structure with all the results copied into it plus their count at the begining
      t_result * rv;
      unsigned int results_count=results.size();
      rv=(t_result *) malloc(results_count*sizeof(t_result));
      memcpy(rv,results.data(),results_count*sizeof(t_result));
      //clean results
      results.clear();
      //unlock results
      results_mutex.unlock();
      //write results into the socket
      write(sockfd,&results_count,sizeof(results_count));
      write(sockfd,rv,results_count*sizeof(t_result));
      free(rv);
      //PHASE 3: GET TARGETS
      //read target count (if -1, set processing_finished)
      int target_count;
      t_pair * tv;
      read(sockfd,&target_count,sizeof(target_count));
      if (target_count==-1){
        processing_finished=true;
      }
      //read all targets in a single structure in memory
      if (target_count>0){
        tv=(t_pair *)malloc(target_count*sizeof(t_pair));
        read(sockfd,tv,target_count*sizeof(t_pair));
        close(sockfd);
        //lock targets
        targets_mutex.lock();
        //insert targets
        for(int i=0;i<target_count;i++) targets.push(tv[i]);
        //unlock targets
        targets_mutex.unlock();
        free(tv);
      } else {
        close(sockfd);
      } 
      //if (target_count>=0) cout<<target_count<<" targets added"<<endl;
      total_targets+=target_count;
      if (debug_max_targets && total_targets>=debug_max_targets) processing_finished=true;
    }
    else {
      targets_mutex.unlock();
      sleep(10);
    }
  }
  cout<< "== Processing finished, waiting for the slaves to die =="<<endl;
  cout<<"TIME SPENT WORKING: "<<(time(NULL) - time_start)<<endl;
  for (int wt=0;wt<threads;wt++) {
    workers[wt].join();
  }

}
