#ifndef FORCE_DEBUG
#define NDEBUG
#endif


#include "WorkQueue.h"


WorkQueue::WorkQueue(int _minLen, string _sQuery, int _queryChunk, string _sTarget, int _targetChunk, double _minProb, double _sigCutoff, bool _probtable, int _slave_count, int _threads){
  minLen=_minLen;
  queryChunk=_queryChunk;
  query_filename=_sQuery;
  targetChunk=_targetChunk;
  target_filename=_sTarget;
  minProb=_minProb;
  probTable=_probtable;
  sigCutoff=_sigCutoff;
  slave_count=_slave_count;
  threads=_threads;
  gethostname(master_hostname, sizeof(master_hostname));
  shutdown_status=0;
  port=MYPORT;
  last_connections.resize(_slave_count);
  slaves_finished_count=0;
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
  pairs_mutex.lock();
  pairs.push_back(p);
  pairs_mutex.unlock();
}
unsigned int WorkQueue::pending_pair_count(){
  unsigned int p=0;
  pairs_mutex.lock();
  for (unsigned int i=0;i<pairs.size();i++)
    if (pairs[i].status==0) p++;
  pairs_mutex.unlock();
  return p;
}

int WorkQueue::next_pair_to_process(unsigned int slave_id){
  int p=0;
  pairs_mutex.lock();
  //TODO: optimization, save the number of the last completely processed pair and start looking from there onwards.
  for (p=0;p<pairs.size()&&pairs[p].status!=0;p++);
  if (p<pairs.size()) {
    pairs[p].slave_id=slave_id;
    pairs[p].status=1;
  } else p=-1;
  pairs_mutex.unlock();
  return p;
}

void WorkQueue::get_pair(t_pair * p, unsigned int pair_id){
  pairs_mutex.lock();
  memcpy(p,&pairs[pair_id],sizeof(t_pair));
  pairs_mutex.unlock();
}

//Accepts connection then returns new socket
//reads client id and command
//TODO: update client_id_last_connection
int WorkQueue::accept_client() {
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
      //touch the slave last connection
      last_connections[conn_clientID-1]==time(NULL);
      return 0;
    }
  }
  close(conn_sockfd);
  conn_sockfd=0;
  conn_clientID=-1;
  return -1;
}

void WorkQueue::submit_tasks() {
  //ALG: caculate how many pairs to assign
  //XXX:number of pairs per slave is fixed
  int tp=0;
  int np=0;
  //XXX: i keep using malloc... can't help it!
  t_pair *p;
  p=(t_pair *)malloc(threads*2*sizeof(t_pair));
  while (tp<threads*2){
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
}

//TODO: implement some kind of check on the reeived data
// if any of the solutions do not check, discard all and mark jobs assigned to this client to be redone
void WorkQueue::receive_solutions() {

  //ALG: reads ammount of results
  //TODO:should check for ammounts of bytes readed here too!!!
  //cout<<"Slave is sending solutions..."<<endl;
  results_mutex.lock();
  slaves_finished_count++;
  results_mutex.unlock();
  unsigned int solcount;
  read(conn_sockfd,&solcount,sizeof(solcount));
  //ALG: reads results into temp variable TODO: validation!!!
  t_result r;
  //cout<<" Current size of resultsv is "<<resultsv.size()<<endl;
  //TODO: read results in a temporary vector<result_t> and only mix aftger all are ok, otherwise we need to recalculate
  for (unsigned int i=0;i<solcount;i++) {
    int availb;
    ioctl(conn_sockfd, FIONREAD, &availb);
    for(int tries=0;tries<500 && availb<sizeof(t_result);tries++){
      //cout<<"Waiting for the slow slave "<< conn_clientID <<"to send the results..."<<endl;
      usleep(10000);
      ioctl(conn_sockfd, FIONREAD, &availb);
    }
    int n=read(conn_sockfd,&r,sizeof(t_result));
    results_mutex.lock();
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
    results_mutex.unlock();
  }
}

void WorkQueue::serve(){
  //socklen_t clilen;
  //struct sockaddr_in serv_addr, cli_addr;
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
      accept_client(); 
      if (conn_sockfd>0){
        receive_solutions();
        submit_tasks();
      }
      close(conn_sockfd);
    }
    //TODO: housekeeping
    //check for stalled slaves and reassign their work.
    //check for the socket to be well and accepting connections.
  }
  //ALG: shutdown
  cout<<"Shutting down all slaves";
  int alive_slaves=slave_count;
  while (alive_slaves){

    //ALG: while slaves alive
    //ALG: wait on socket
    struct timeval tv;
    tv.tv_sec=1;
    tv.tv_usec=0;
    FD_ZERO(&rfds);
    FD_SET(sockfd, &rfds);
    if (select(sockfd+1, &rfds, (fd_set *) 0, (fd_set *) 0, &tv)>0){
      accept_client(); 
      if (conn_sockfd>0){
        receive_solutions();
        //ALG: I don't care about you, just DIE
        int tp=-1;
        write(conn_sockfd,&tp,sizeof(tp));
        alive_slaves--;
      }
      close(conn_sockfd);
    }
  }
  return;
  //ALG: ---Listen Function Ends
}

void WorkQueue::join(){
  server_thread.join();
}


void WorkQueue::start_listener(){
  //TODO opens socket
  int portno;
  //socklen_t clilen;
  //char buffer[256];
  struct sockaddr_in serv_addr;//, cli_addr;
  int n;
  sockfd = socket(AF_INET, SOCK_STREAM, 0);

  if (sockfd < 0) {
    cout<<"ERROR opening socket"<<endl;
    exit(1);
  }
  //SET TO BLOCKING!!!!!
  //int flags = fcntl(sockfd, F_GETFL);
  //int result = fcntl(sockfd, F_SETFL, flags & ~O_NONBLOCK);
  memset((char *) &serv_addr, '\0', sizeof(serv_addr));
  portno = port;
  serv_addr.sin_family = AF_INET;
  serv_addr.sin_addr.s_addr = INADDR_ANY;
  //try up to 100 port numbers (in case port is used)
  for(portno=port;portno<port+100;portno++){
    serv_addr.sin_port = htons(portno);
    if (::bind(sockfd, (struct sockaddr *) &serv_addr, sizeof(serv_addr)) < 0) 
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
  server_thread=std::thread(&WorkQueue::serve,this);
}

void WorkQueue::close_queue(){
  shutdown_status=1;
}

void WorkQueue::results_from_file(const char * filename, const std::vector<SeqChunk> & _queryInfo){
  FILE * rf=fopen(filename,"r");
  //preallocate vector space
  fseek(rf, 0 , SEEK_END);
  std::cout<<"Pre-Allocating memory for "<<ftell(rf)/sizeof(t_result)<<" new elements"<<std::endl;
  resultsv.reserve(resultsv.size()+ftell(rf)/sizeof(t_result));
  std::cout<<"Done!"<<std::endl;
  fseek(rf, 0 , SEEK_SET);
  t_result r;
  while (fread(&r,sizeof(r),1,rf)==1){
    if (r.reverse) {//recalculate postion in satsuma format
      //r.qstart=r.qstart + r.query_size - _queryInfo[r.query_id].GetStart() - queryChunk;
      r.qstart=r.query_size-r.len-r.qstart;
    }
    resultsv.push_back(r);
  }
  fclose(rf);
}

void WorkQueue::setup_queue(const int mem){
  //TODO: avoid all hardcoding and support PBS/LSF
  //spawns each slave with its slave_id
  start_listener();  
  stringstream cmd;
  stringstream sh_cmd;
  for (int i=0;i<slave_count;i++){
    cmd.str("");
    cmd << std::getenv("SATSUMA2_PATH") << "/HomologyByXCorrSlave" << " -master " << master_hostname << " -port " << port;
    cmd << " -sid " << i+1 << " -p "<< threads << " -q " << query_filename << " -t " << target_filename;
    cmd << " -l " << minLen << " -q_chunk " << queryChunk << " -t_chunk " << targetChunk << " -min_prob " << minProb;
    cmd << " -cutoff " << sigCutoff << (probTable ? " -prob_table true": "");

    cout<< "Launching slave: " << endl << "  " << cmd.str() <<endl;

    // JW: slaves will always be run asynchronously
    sh_cmd.str("");
    sh_cmd << "sh " << std::getenv("SATSUMA2_PATH") << "/satsuma_run.sh " << std::getenv("PWD") << " \"" << cmd.str() << "\" " << to_string(threads);
    sh_cmd << " " << to_string(mem) << " SL" << i+1 << " " << to_string(0);

    system(sh_cmd.str().c_str());
  }
}


t_collect_status WorkQueue::collect_new_matches(MultiMatches &matches){
  //processes new matches and returns hom many there where
  //ALG: mutex start
  unsigned long int i;
  t_collect_status status;
  unsigned long int revcount=0,fwcount=0;
  results_mutex.lock();
  matches.reserve(matches.GetMatchCount()+resultsv.size());
  //ALG: update the MultiMatch
  SingleMatch m;
  for (i=0;i<resultsv.size();i++){
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
    if (resultsv[i].reverse) revcount++;
    else fwcount++;
  }
  //ALG: mutex end
  resultsv.clear();
  status.slaves=slaves_finished_count;
  status.matches=i;
  if (status.matches>0 and status.slaves==0) status.slaves=1; //XXX: wrong way to resolve race!!!
  slaves_finished_count=0;
  results_mutex.unlock();
  cout<<"WORKQUEUE: matches collected. FW: "<<fwcount<<"   REV: "<<revcount<<endl;
  return status;
}
