#ifndef FORCE_DEBUG
#define NDEBUG
#endif


#include <string.h>
#include <unistd.h>
#include <fcntl.h> 
#include <ctime>
#include "base/CommandLineParser.h"
#include "analysis/SequenceMatch.h"
#include "analysis/GridSearch.h"
#include "analysis/SeqChunk.h"
#include "analysis/MatchDynProg.h"
#include "util/SysTime.h"
#include <math.h>
#include <pthread.h>
#include <netinet/in.h> 
#include <sys/socket.h>
#include <sys/ioctl.h>
#include "analysis/WorkQueue.h"


WorkQueue::WorkQueue(int _minLen, string _sQuery, int _queryChunk, string _sTarget, int _targetChunk, double _minProb, double _sigCutoff, int _slave_count){
  minLen=_minLen;
  queryChunk=_queryChunk;
  query_filename=_sQuery;
  targetChunk=_targetChunk;
  target_filename=_sTarget;
  minProb=_minProb;
  sigCutoff=_sigCutoff;
  slave_count=_slave_count;
  gethostname(master_hostname, sizeof(master_hostname));
  shutdown_status=0;
  port=MYPORT;
  pthread_mutex_init(&pairs_mutex,NULL);
  last_connections.resize(_slave_count);
}

void WorkQueue::add_pair(int _targetFrom, int _targetTo, int _queryFrom, int _queryTo, bool _fast){
  t_pair p;
  p.targetFrom=_targetFrom;
  p.targetTo=_targetTo;
  p.queryFrom=_queryFrom;
  p.queryTo=_queryTo;
  p.fast=_fast;
  p.slave_id=0;
  p.status=0;
  pthread_mutex_lock(&pairs_mutex);
  pairs.push_back(p);
  pthread_mutex_unlock(&pairs_mutex);
}
unsigned int WorkQueue::pending_pair_count(){
  unsigned int p=0;
  pthread_mutex_lock(&pairs_mutex);
  for (unsigned int i=0;i<pairs.size();i++)
    if (pairs[i].status==0) p++;
  pthread_mutex_unlock(&pairs_mutex);
  return p;
}

int WorkQueue::next_pair_to_process(unsigned int slave_id){
  int p=0;
  pthread_mutex_lock(&pairs_mutex);
  //TODO: optimization, save the number of the last completely processed pair and start looking from there onwards.
  for (p=0;p<pairs.size()&&pairs[p].status!=0;p++);
  if (p<pairs.size()) {
    pairs[p].slave_id=slave_id;
    pairs[p].status=1;
  } else p=-1;
  pthread_mutex_unlock(&pairs_mutex);
  return p;
}

void WorkQueue::get_pair(t_pair * p, unsigned int pair_id){
  pthread_mutex_lock(&pairs_mutex);
  memcpy(p,&pairs[pair_id],sizeof(t_pair));
  pthread_mutex_unlock(&pairs_mutex);
}

static void * thread_serve(void * args){
  //ALG: ---Listen Function Starts
  WorkQueue * wq= (WorkQueue *)args;
  wq->serve();
}

//Accepts connection then returns new socket
//reads client id and command
//TODO: update client_id_last_connection
int WorkQueue::accept_client_command() {
  //cout<<"***Incoming connection!!!****"<<endl;
  conn_sockfd = accept(sockfd, NULL, NULL); //just discard the client info!
  if (conn_sockfd < 0) {
    cout<< "ERROR on accept." <<endl;
  } else {
    //ALG: read slave_id
    int n = read(conn_sockfd,&conn_clientID,sizeof(conn_clientID));
    if (n!=sizeof(conn_clientID) || conn_clientID > slave_count){
      cout<<"ERROR No valid slave ID received.";
    } else {
      int reply_ok=STATUS_OK;
      write(conn_sockfd,&reply_ok,sizeof(reply_ok));

      n = read(conn_sockfd,&conn_command,sizeof(conn_command));
      if (n!=sizeof(conn_command) || (conn_command!=TCP_COMMAND_REQUEST_PAIRS && conn_command!=TCP_COMMAND_SEND_SOLUTIONS)){
        cout<<"ERROR: no valid command received from slave "<<conn_clientID<<endl;
      } else {
        //touch the slave last connection
        last_connections[conn_clientID-1]==time(NULL);
        return 0;
      }
    }
  }
  close(conn_sockfd);
  conn_sockfd=0;
  conn_clientID=-1;
  conn_command=0;
  return -1;
}

void WorkQueue::submit_tasks() {
  //ALG: caculate how many pairs to assign
  //XXX:number of pairs per slave is fixed
  int tp=0;
  int np=0;
  //XXX: i keep using malloc... can't help it!
  t_pair *p;
  p=(t_pair *)malloc(4*sizeof(t_pair));
  while (tp<4){
    np=next_pair_to_process(conn_clientID);
    //cout<<"  next pair to process is number "<<np<<endl;
    if (np==-1) break;
    get_pair(&p[tp],np);
    tp++;
  }
  write(conn_sockfd,&tp,sizeof(tp));
  write(conn_sockfd,p,sizeof(t_pair)*tp);
  if (tp){
    //cout<<"Slave "<<conn_clientID<<" will process "<<tp<<" pairs, I have "<<pending_pair_count()<<" pending pairs on the queue"<<endl;
    if (idle_slaves[conn_clientID]) {
      cout<<"Slave "<<conn_clientID<<" is working now."<<endl;
      idle_slaves[conn_clientID]=false;
    }
  } else {
    if (!idle_slaves[conn_clientID-1]) {
      idle_slaves[conn_clientID-1]=true;
      cout<<"Slave "<<conn_clientID<<" will remain idle."<<endl;
    }
  }
  free(p);
  close(conn_sockfd);
}

//TODO: implement some kind of check on the reeived data
// if any of the solutions do not check, discard all and mark jobs assigned to this client to be redone
void WorkQueue::receive_solutions() {

  //ALG: reads ammount of results
  //TODO:should check for ammounts of bytes readed here too!!!
  //cout<<"Slave is sending solutions..."<<endl;
  unsigned int solcount;
  read(conn_sockfd,&solcount,sizeof(solcount));
  //ALG: reads results into temp variable TODO: validation!!!
  t_result r;
  //cout<<" Current size of resultsv is "<<resultsv.size()<<endl;
  //TODO: read results in a temporary vector<result_t> and only mix aftger all are ok, otherwise we need to recalculate
  for (unsigned int i=0;i<solcount;i++) {
    int availb;
    ioctl(conn_sockfd, FIONREAD, &availb);
    for(int tries=0;tries<50 && availb<sizeof(t_result);tries++){
      //cout<<"Waiting for the slow slave to send the results..."<<endl;
      usleep(100000);
      ioctl(conn_sockfd, FIONREAD, &availb);
    }
    int n=read(conn_sockfd,&r,sizeof(t_result));
    pthread_mutex_lock(&results_mutex);
    resultsv.push_back(r);
    if (r.len==0 || r.query_size==0 || n!=sizeof(t_result)) {cout<<"ERROR: NEW Match received from slave "<<conn_clientID<<" saved to position "<<resultsv.size()<<" looks wrong!"<<endl;
      cout<<": "<<n<<"/"<<sizeof(t_result)<<" bytes received"<<endl;
      cout<<": query_id="<<r.query_id<<endl;
      cout<<": target_id="<<r.target_id<<endl;
      cout<<": query_size="<<r.query_size<<endl;
      cout<<": qstart="<<r.qstart<<endl;
      cout<<": tstart="<<r.tstart<<endl;
      cout<<": len="<<r.len<<endl;
      cout<<": reverse="<<r.reverse<<endl;
      cout<<": prob="<<r.prob<<endl;
      cout<<": ident="<<r.ident<<endl;
    }
    pthread_mutex_unlock(&results_mutex);
  }
  int reply_ok=STATUS_OK;
  write(conn_sockfd,&reply_ok,sizeof(reply_ok));
  close(conn_sockfd);

}

void WorkQueue::serve(){
  socklen_t clilen;
  struct sockaddr_in serv_addr, cli_addr;
  int n;
  for (int i=0;i<slave_count;i++) idle_slaves.push_back(false);
  //ALG: while the shutdown flag is not true
  listen(sockfd,1000);//TODO: accept more than 1000 requests? how many children are we having?
  fd_set rfds;

  while (!shutdown_status){
    //ALG: select on socket (timeout of a second to allow shutdown?)
    struct timeval tv;
    tv.tv_sec=1;
    tv.tv_usec=0;
    FD_ZERO(&rfds);
    FD_SET(sockfd, &rfds);
    if (select(sockfd+1, &rfds, (fd_set *) 0, (fd_set *) 0, &tv)>0){
      accept_client_command(); //If connection is closed, the command will be 0, so nothing more happens
      if (conn_command==TCP_COMMAND_REQUEST_PAIRS){
        submit_tasks();
      }
      else if (conn_command==TCP_COMMAND_SEND_SOLUTIONS){
        receive_solutions();
      }
    }
    //TODO: housekeeping
    //check for stalled slaves and reassign their work.
    //check for the socket to be well and accepting connections.
  }
  //ALG: shutdown
  cout<<"Shutting down all slaves";
  int alive_slaves=slave_count;
  while (slave_count){

    //ALG: while slaves alive
    //ALG: wait on socket
    //ALG: I don't care about you, just DIE
  }
  return;
  //ALG: ---Listen Function Ends
}


void WorkQueue::start_listener(){
  //TODO opens socket
  int portno;
  socklen_t clilen;
  char buffer[256];
  struct sockaddr_in serv_addr, cli_addr;
  int n;
  sockfd = socket(AF_INET, SOCK_STREAM, 0);

  if (sockfd < 0) {
    cout<<"ERROR opening socket"<<endl;
    exit(1);
  }
  //SET TO BLOCKING!!!!!
  int flags = fcntl(sockfd, F_GETFL);
  int result = fcntl(sockfd, F_SETFL, flags & ~O_NONBLOCK);
  memset((char *) &serv_addr, '\0', sizeof(serv_addr));
  portno = port;
  serv_addr.sin_family = AF_INET;
  serv_addr.sin_addr.s_addr = INADDR_ANY;
  //try up to 100 port numbers (in case port is used)
  for(portno=port;portno<port+100;portno++){
    serv_addr.sin_port = htons(portno);
    if (bind(sockfd, (struct sockaddr *) &serv_addr,
          sizeof(serv_addr)) < 0) 
      cout<<"ERROR on binding"<<endl;
    else
      break;
  }
  if (portno==port+100){ //No binding was done!
    cout<<"Run out of ports... sorry!"<<endl;
    exit(1);
  }
  //TODO spawns a thread to listen on the socket
  port=portno;
  cout<<"Spawning a thread to serve the WorkQueue on port "<<portno<<endl;
  pthread_create(&server_thread,NULL,thread_serve,this);
}

void WorkQueue::close_queue(){
  shutdown_status=1;
}

void WorkQueue::setup_queue(){
  //TODO: avoid all hardcoding and support PBS/LSF
  //spawns each slave with its slave_id
  start_listener();  
  stringstream cmd;
  cmd << "echo '" << std::getenv("SATSUMA2_PATH") << "/HomologyByXCorrSlave";
  for (int i=0;i<slave_count;i++){
    cmd << " -master " << master_hostname << " -port " << port << " -sid " << i+1;
    cmd << " -q " << query_filename << " -t " << target_filename;
    cmd << " -l " << minLen << " -q_chunk " << queryChunk << " -t_chunk " << targetChunk << " -min_prob " << minProb << " -cutoff " << sigCutoff << " &";
    if (i%8==7 || i==slave_count-1){
      cmd << " wait '|qsub -l ncpus=" << 8 << " -N SL" << i+1;
      cout<< "Launching slave with command line:"<<endl<<"  "<<cmd.str()<<endl;
      system(cmd.str().c_str());
      cmd.str("");
      cmd << "echo '" << std::getenv("SATSUMA2_PATH") << "/HomologyByXCorrSlave";
    }
  }
}


unsigned int WorkQueue::collect_new_matches(MultiMatches &matches){
  //processes new matches and returns hom many there where
  //ALG: mutex start
  unsigned long int i;
  pthread_mutex_lock(&results_mutex);
  //ALG: update the MultiMatch
  for (i=0;i<resultsv.size();i++){
    SingleMatch m;
    m.SetQueryTargetID(resultsv[i].query_id, resultsv[i].target_id, resultsv[i].query_size);
    m.SetPos(resultsv[i].qstart, resultsv[i].tstart, resultsv[i].len, resultsv[i].reverse);
    m.SetProbability(resultsv[i].prob);
    m.SetIdentity(resultsv[i].ident);
    m.AddMatches(resultsv[i].ident * (double)resultsv[i].len);
    if (m.GetLength()==0) {
      cout<<"ERROR!!! Match with Length==0 collected (resultsv["<<i<<"])"<<endl;
      cout<<"query_id="<<resultsv[i].query_id<<endl;
      cout<<"target_id="<<resultsv[i].target_id<<endl;
      cout<<"query_size="<<resultsv[i].query_size<<endl;
      cout<<"qstart="<<resultsv[i].qstart<<endl;
      cout<<"tstart="<<resultsv[i].tstart<<endl;
      cout<<"len="<<resultsv[i].len<<endl;
      cout<<"reverse="<<resultsv[i].reverse<<endl;
      cout<<"prob="<<resultsv[i].prob<<endl;
      cout<<"ident="<<resultsv[i].ident<<endl;
    }
    matches.AddMatch(m);
  }
  //ALG: mutex end
  resultsv.clear();
  pthread_mutex_unlock(&results_mutex);
  return i;
}
