#ifndef _WORKQUEUE_H_
#define _WORKQUEUE_H_

#include <string.h>
#include <unistd.h>
#include <fcntl.h>
#include "analysis/SequenceMatch.h"
#include "analysis/SeqChunk.h"
#include "util/SysTime.h"
#include <math.h>
#include <thread>
#include <mutex>
#include <netinet/in.h>
#include <sys/socket.h>
#include <sys/ioctl.h>

typedef struct {
  int targetFrom, targetTo, queryFrom, queryTo;
  bool fast;
  int slave_id;
  unsigned char status; 
} t_pair;
typedef struct {
    unsigned long int query_id;
    unsigned long int target_id;
    unsigned long int query_size;
    unsigned long int qstart;
    unsigned long int tstart;
    unsigned long int len;
    bool reverse;
    double prob;
    double ident;
} t_result;
typedef struct {
  uint64_t slaves;
  uint64_t matches;
} t_collect_status;

//XXX: this belongs in a lib
#define MYPORT 3491
#define STATUS_OK 100
#define TCP_COMMAND_REQUEST_PAIRS 200
#define TCP_COMMAND_SEND_SOLUTIONS 202
#define TIMEOUT_SEC 3

#define PAIR_PENDING 0
#define PAIR_ASSIGNED 1
#define PAIR_COMPUTED 2
#define PAIR_DONE 3


//TODO: make the serve function a method, use a wrapper to call the method from the thread execution
class WorkQueue {
  public:
    WorkQueue(int _minLen, string _sQuery, int _queryChunk, string _sTarget, int _targetChunk, double _minProb, double _sigCutoff, int _slave_count, int _threads);
    void results_from_file(const char * ,const svec<SeqChunk> & _queryInfo);
    void setup_queue();
    void add_pair(int targetFrom, int targetTo, int queryFrom, int queryTo, bool fast);
    unsigned int pending_pair_count();
    t_collect_status collect_new_matches(MultiMatches &matches); //updates the multimatches and returns number of new matches
    void close_queue();
    void serve();//spawns a new thread and serves the Queue

  private:
    void start_listener();
    int next_pair_to_process(unsigned int slave_id);//gets a pair to process, and marks it
    void get_pair(t_pair * p, unsigned int pair_id);

    int accept_client();
    void submit_tasks();
    void receive_solutions();

    char shutdown_status; // 0-normal operation; 1-shutdown requested; 2-clients finished and thread exited;
    int slave_count;
    int threads;
    std::thread server_thread;
    int minLen, queryChunk, targetChunk;
    string query_filename,target_filename;
    double minProb, sigCutoff;
    char master_hostname[1024];//XXX: hardcoded!!!
    int port;
    unsigned long int slaves_finished_count;

    std::mutex pairs_mutex;
    std::vector<t_pair> pairs; //XXX: beware, all access to this should be mutexed!!!
    std::mutex results_mutex;
    std::vector<t_result> resultsv;
    
    int sockfd;
    int conn_sockfd;//socket for the connection being processed XXX: this limits to 1 slave being served at a time.
    int conn_clientID;
    int conn_command;
    std::vector<time_t> last_connections;
    std::vector<bool> idle_slaves;
    

};

#endif
