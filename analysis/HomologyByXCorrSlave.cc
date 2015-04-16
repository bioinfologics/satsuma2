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
#include "util/mutil.h"
//#include "base/FileParser.h"
#include <pthread.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <netdb.h>


//====== XXX: REALLY????? =======================================
int RCQuery(int offset, int start, int len, int chunkSize, int qLen) 
{
  int l = start + qLen - offset - chunkSize;
  return l;
}

//void Load(vecDNAVector & out, svec<string> & names, const string & file) 
//{
//  out.Read(file, names);

//}
//==============================================================

class HomologyByXCorr
{
public:
  HomologyByXCorr(string mHostname, int port, int _sid, string qFilename, int qChunk, string tFilename, int tChunk, double _topCutoff, double _topCutoffFast) {
    m_minLen = 0;
    m_pMulti = NULL;
    m_offset = 0;
    m_qLen = 0;
    m_chunk = 0;
    m_tOffset = 0;
    m_qID = -1;
    m_tID = -1;

    m_targetSize = 0;
    m_minProb = 0.99;
    queryChunk=qChunk;
    targetChunk=tChunk;
    query_filename=qFilename;
    target_filename=tFilename;
    master_hostname=mHostname;
    portno=port;
    slave_id=_sid;
    pairs=NULL;
    topCutoff=_topCutoff;
    topCutoffFast=_topCutoffFast;
  }

  void SetMinProb(double p) {
    m_minProb = p;
  }

  void Align(int _target_id, CCSignal _target_signal, int _query_id, CCSignal _query_signal, CCSignal _query_rcsignal, bool fast);

  void SetTargetSize(double i) {
    m_targetSize = i;
  }

  void SetMinimumAlignLen(int l) {
    m_minLen = l;
  }

  void disconnect_from_master();
  int connect_to_master();
  void create_chunks();
  int get_targets_from_master();
  void send_solutions_to_master();
  std::vector<CCSignal> create_signals(vecDNAVector _sequences, unsigned long int _from, unsigned long int _count, bool reverse);
  void align_block(unsigned long int pair_id);
  void FilterMatches(int _target_id, int _query_id, vecSeqMatch _matches, bool _reverse);
  void work();


private:
  int sid;
  t_pair * pairs;
  int m_minLen;
  double topCutoff, topCutoffFast;
  double targetTotal;
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
  int queryChunk,targetChunk;
  string target_filename,query_filename;
  vecDNAVector target, query;
  svec<string> targetNames, queryNames; //holds the sequences names
  ChunkManager * cmQuery; //XXX: forced to use this as no default init
  ChunkManager * cmTarget; //XXX: forced to use this as no default init
  svec<SeqChunk> targetInfo, queryInfo;
  std::vector<t_result> resultsv;
 


};

void HomologyByXCorr::create_chunks(){
  //Loads both the query and the target files into memory
  vecDNAVector targetRaw, queryRaw; //holds the sequences
  cout << "Loading query sequence:  " << query_filename << endl;
  queryRaw.Read(query_filename,queryNames);
  cout << " - Creating chunks..."<<endl;
  cmQuery = new ChunkManager(queryChunk,0);
  cmQuery->ChunkIt(query, queryInfo, queryRaw, queryNames, 0, 0);
  queryRaw.clear();


  cout << "Loading target sequence:  " << target_filename << endl;
  targetRaw.Read(target_filename,targetNames);
  cout << " - Creating chunks..."<<endl;
  cmTarget = new ChunkManager(targetChunk, targetChunk / 4);
  cmTarget->ChunkIt(target, targetInfo, targetRaw, targetNames, 0, 0);
  targetRaw.clear();
  cout << "Done loading and chunking." << endl;

  //TODO: set total target size and check if needed to create the multi!!!!!!
  targetTotal = 0;
  for (int i=0; i<cmTarget->GetCount(); i++) {
    targetTotal += (double)cmTarget->GetSize(i);
  }

}


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
  
  if (write(sockfd,&slave_id,sizeof(slave_id)) != sizeof(slave_id)) {
    cout << "ERROR sending id to master" << endl;
    close(sockfd);
    return -1;
  }
  int server_reply;
  if (read(sockfd,&server_reply,sizeof(server_reply)) < 0) {
    cout << "ERROR reading reply from master" << endl;
    close(sockfd);
    return -1;
  }
  if (server_reply!=STATUS_OK){
    cout << "ERROR master did not reply OK" << endl;
    close(sockfd);
    return -1;
  }
  
  return sockfd;
    
}

int HomologyByXCorr::get_targets_from_master(){
	//connects to master, gets targets, if master says no targets, wait TIMEOUT_SEC seconds and try again
	//Sets the targets array, and the num_targets counter
  //If master say DIE: num_targets=0
  // Create socket, if can't connect
  int pair_count;
  while(1){
    sockfd=connect_to_master();
    if (sockfd<=0) {
      cout<<"ERROR: can't connect to master! will retry later..."<<endl;
      //return -1;
      sleep(5);
      continue;
    }
    //Send request for new PAIRS
    //TODO:write(TCP_COMMAND_REQUEST_PAIRS);
    //TODO: read pair count
    unsigned int command;
    command=TCP_COMMAND_REQUEST_PAIRS;
    write(sockfd,&command,sizeof(command));
    if (read(sockfd,&pair_count,sizeof(pair_count)) < sizeof(pair_count)){
      cout<<"ERROR: failure on reading socket"<<endl;
      return -1;
    }

    if (pair_count>0){//Master has work for us, let's read it
      cout<<"Receiving "<<pair_count<<" pairs from master"<<endl;
      //TODO: malloc the pairs array
      if (pairs!=NULL) free(pairs);
      pairs=(t_pair *)malloc(pair_count*sizeof(t_pair));
      //TODO: read pairs
      read(sockfd,pairs,pair_count*sizeof(t_pair));
    }
    disconnect_from_master();
    if (pair_count==-1){
      return -1; //master said DIE
    }
    if (pair_count==0){//Master has not any work, disconnect and wait
      //TODO: wait
      sleep(TIMEOUT_SEC);
    }
    if (pair_count>0) return pair_count;
  }
}

void HomologyByXCorr::send_solutions_to_master(){
    //connects to master, sends solutions and clears them
    int sockfd;
    int pair_count;
    unsigned int command;
    while(1){
        sockfd=connect_to_master();
        if (sockfd<=0) {
            cout<<"ERROR: can't connect to master! (Retrying in 1 second)"<<endl;
            sleep(1);
        }
        else break;
    }
    cout<<" Sending "<<resultsv.size()<<" solutions"<<endl;
    //Send request to upload solutions
    command=TCP_COMMAND_SEND_SOLUTIONS;
    write(sockfd,&command,sizeof(command));
    //write solution count
    unsigned int solcount=resultsv.size();
    write(sockfd,&solcount,sizeof(solcount));
    //write solutions
    for (unsigned int i=0;i<solcount;i++) write(sockfd,&resultsv[i],sizeof(t_result));
    int reply_ok=0;
    read(sockfd,&reply_ok,sizeof(reply_ok));
    if (reply_ok==STATUS_OK){
        cout<<" Server finished reading OK."<<endl;
    }
    disconnect_from_master();
    //TODO: free(pairs);
    free(pairs);
    resultsv.resize(0);
    pairs=NULL;

    pair_count=0;
    cout<<" SENT"<<endl;
}

void HomologyByXCorr::FilterMatches(int _target_id, int _query_id, vecSeqMatch _matches, bool _reverse){
    int len, tStart, qStart;
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

        double prob = GetMatchProbability(target[_target_id],
                query[_query_id],
                _matches[j].GetStartTarget(),
                _matches[j].GetStartQuery(),
                _matches[j].GetLength(),
                targetTotal);

        if (prob < m_minProb)
            continue;

        //cout << "Match # " << j << " probability " << 100.* prob << " %" << endl;
        //cout << "Start target: " << tStart << " - " << tStart + len << endl;
        double ident = PrintMatch(query[_query_id], target[_target_id], _matches[j], true);//XXX: does this really print? PERFORMANCE!!!
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
        resultsv.push_back(r);
    }

}

void HomologyByXCorr::Align(
        int _target_id,
        CCSignal _target_signal, 
        int _query_id,
        CCSignal _query_signal, 
        CCSignal _query_rcsignal, 
        bool _fast)
{
    CrossCorrelation xc;
    vecSeqMatch matches,revmatches;
    SeqAnalyzer sa;
    int i, j;
    svec<float> result,revresult; 

    sa.SetTopCutoff((_fast ? topCutoffFast : topCutoff));

    xc.CrossCorrelate(result, _target_signal, _query_signal);
    matches.clear();
    sa.MatchUp(matches, query[_query_id], target[_target_id], result);
    FilterMatches(_target_id, _query_id, matches, false);

    xc.CrossCorrelate(revresult, _target_signal, _query_rcsignal);
    revmatches.clear();
    sa.MatchUp(revmatches, query[_query_id], target[_target_id], revresult);
    FilterMatches(_target_id, _query_id, revmatches, true);

}
std::vector<CCSignal> HomologyByXCorr::create_signals( vecDNAVector _sequences, unsigned long int _from, unsigned long int _count, bool reverse){
    cout<<"creating signals..."<<endl;
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

void HomologyByXCorr::align_block(unsigned long int pair_id){

    cout<<"Align block running (block: "<<pair_id<<" - t:"<<pairs[pair_id].targetFrom<<"-"<<pairs[pair_id].targetTo<<" q:"<<pairs[pair_id].queryFrom<<"-"<<pairs[pair_id].queryTo<<")"<<endl;
    //XXX: is this storage correct? create all signals for the target and query spaces
    std::vector<CCSignal> targetSignals;
    std::vector<CCSignal> querySignals;
    std::vector<CCSignal> querySignalsReverse;
    unsigned long int target_count, query_count;

    //Create all signals
    //XXX: is To included or not?
    target_count=pairs[pair_id].targetTo+1-pairs[pair_id].targetFrom;
    query_count=pairs[pair_id].queryTo+1-pairs[pair_id].queryFrom;
    //XXX what is the coordinates, how are them set?
    targetSignals = create_signals(target, pairs[pair_id].targetFrom, target_count, false);
    querySignals = create_signals(query, pairs[pair_id].queryFrom, query_count, false);
    querySignalsReverse = create_signals(query, pairs[pair_id].queryFrom, query_count, true);

    cout<<"Block "<<pair_id<<" signals created..."<<endl;

    for (unsigned long int qi=0; qi<query_count; qi++){
        //for each query
        for (unsigned long int ti=0; ti<target_count; ti++){
            //for each target
            Align(
                    pairs[pair_id].targetFrom+ti, targetSignals[ti],
                    pairs[pair_id].queryFrom+qi, querySignals[qi], querySignalsReverse[qi],
                    pairs[pair_id].fast);
        }
    }
}

void HomologyByXCorr::work(){
    //Work loop, gets targets from master and reports back
    while(1){
        //======= Connect to master, get instructions =======
        cout<<"Getting targets from master"<<endl;
        num_targets=get_targets_from_master();

        cout<<"DONE, got "<<num_targets<<" targets"<<endl;
        if (0>=num_targets) break; //Master said DIE!!!

        cout<<"Aligning..."<<endl;
        for (int i=0;i<num_targets;i++) align_block(i);
        //======= send reply to master
        cout<<"Sending solutions"<<endl;
        send_solutions_to_master();
    }
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
    commandArg<string> aStringCmmd("-q","query fasta sequence");
    commandArg<string> bStringCmmd("-t","target fasta sequence");
    commandArg<int> lIntCmmd("-l","minimum alignment length", 0);
    commandArg<int> qChunkCmmd("-q_chunk","query chunk size", 4096);
    commandArg<int> tChunkCmmd("-t_chunk","target chunk size", 4096);
    commandArg<double> probCmmd("-min_prob","minimum probability to keep match", 0.9999);
    commandArg<double> cutoffCmmd("-cutoff","signal selection cutoff", 1.8);
    commandArg<double> cutoffFastCmmd("-cutoff_fast","signal selection cutoff (fast)", 2.9);


    commandLineParser P(argc,argv);
    P.SetDescription("Compares two sequences via cross-correlation");
    P.registerArg(mStringCmmd);
    P.registerArg(portIntCmmd);
    P.registerArg(sidIntCmmd);
    P.registerArg(aStringCmmd);
    P.registerArg(bStringCmmd);
    P.registerArg(lIntCmmd);

    P.registerArg(qChunkCmmd);
    P.registerArg(tChunkCmmd);
    P.registerArg(probCmmd);
    P.registerArg(cutoffCmmd);
    P.registerArg(cutoffFastCmmd);

    P.parse();

    string master = P.GetStringValueFor(mStringCmmd);
    int port = P.GetIntValueFor(portIntCmmd);
    int sid = P.GetIntValueFor(sidIntCmmd);
    string sQuery = P.GetStringValueFor(aStringCmmd);
    string sTarget = P.GetStringValueFor(bStringCmmd);
    int minLen = P.GetIntValueFor(lIntCmmd);
    int targetChunk = P.GetIntValueFor(tChunkCmmd);
    int queryChunk = P.GetIntValueFor(qChunkCmmd);
    double minProb = P.GetDoubleValueFor(probCmmd);

    double topCutoff = P.GetDoubleValueFor(cutoffCmmd);
    double topCutoffFast = P.GetDoubleValueFor(cutoffFastCmmd);

    //======= Initialize Classes =======
    cout << "Initializing HomologyByXCorr class..."<<endl;
    HomologyByXCorr hbxc(master,port,sid,sQuery,queryChunk,sTarget,targetChunk,topCutoff,topCutoffFast);
    hbxc.SetMinimumAlignLen(minLen);
    cout<<"DONE"<<endl;
    //======= Load Genomes =======
    hbxc.create_chunks();

    //======= Main loop
    cout<< "== Going to main work loop==";
    hbxc.work();

}
