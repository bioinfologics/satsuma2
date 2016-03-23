#ifndef THREADHANDLER_H
#define THREADHANDLER_H

using namespace std;

#include <pthread.h>
#include <errno.h>
#include <signal.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <iostream>
#include <vector>
#include <queue>
#include <unistd.h>



class ThreadMutex
{
 public:
  ThreadMutex() {
    int r = pthread_mutex_init(&m_lock, NULL);
    if (r > 0)
      cout << "Thread INIT error: " << r << endl;
  
  }
  ~ThreadMutex() {
    pthread_mutex_destroy(&m_lock);
  }

  void Lock() {
    int r = pthread_mutex_lock(&m_lock);
    if (r > 0)
      cout << "Thread LOCK error: " << r << endl;
  }
  void Unlock() {
    int r = pthread_mutex_unlock(&m_lock);
    if (r > 0)
      cout << "Thread UNLOCK error: " << r << endl;
    
  }

  void Print() const {
    const char * p = (const char*)&m_lock;
    for (int i=0; i<(int)sizeof(m_lock); i++)
      cout <<  (int)p[i] << endl;
  }

 private:
  pthread_mutex_t m_lock;
};


class IOneThread
{
 public:
  IOneThread() {
    m_bDo = false;
    m_bDie = false;
    m_iThreadId = 0;
    m_bFinished = false;
    pthread_mutex_init(&m_lock, NULL);
  }

  virtual ~IOneThread() {
    pthread_mutex_destroy(&m_lock);
  }

  string getName() {
	  return name;
  }

  virtual bool Initialize(const string & msg, const string &name = "");

  bool Die() {
    pthread_mutex_lock(&m_lock);
    m_bDie = true;
    pthread_mutex_unlock(&m_lock);
    return true;
  }

  bool Do(const string & msg) {
    pthread_mutex_lock(&m_lock);
    m_msg.push_back(msg);
    m_bDo = true;
    pthread_mutex_unlock(&m_lock);
    return true;
  }

  bool Done() {
    bool bDone = false;
    pthread_mutex_lock(&m_lock);
    bDone = !m_bDo;
    pthread_mutex_unlock(&m_lock);
    return bDone;
  }

  bool Finished() {
    bool bDone = false;
    pthread_mutex_lock(&m_lock);
    bDone = m_bFinished;
    pthread_mutex_unlock(&m_lock);
    return bDone;
  }

  
  void CallInit();

 protected:
  virtual bool OnDie() = 0;
  virtual bool OnDo(const string & msg) = 0;
  virtual bool OnInitialize(const string & msg) = 0;

  bool ShouldDie() {
    bool bDie = false;
    pthread_mutex_lock(&m_lock);
    bDie = m_bDie;
    pthread_mutex_unlock(&m_lock);
    return bDie;
  }

 private:
  bool m_bDie;
  bool m_bDo;
  bool m_bFinished;
  pthread_mutex_t m_lock;
  vector<string> m_msg;
  string m_init;
  string name;

  pthread_t m_iThreadId;


};


//=================================================

class ThreadHandler
{
 public:
  ThreadHandler();
  ~ThreadHandler();

  // Important: the ThreadHandler will delete the thread objects!!
  int AddThread(IOneThread * pThread, const string & init_msg = "", const string & name = "");

  int GetNumThreads() const {return (int)m_threads.size();}
  int GetThreadID(int i) const {return i;}

  bool Feed(int threadID, const string & msg);
  //Gives the message to the first finished thread
  void FeedAny(const string &msg);

  //Currently, this function has undefined behavior if a thread may be calling FeedAny
  //simultaneously. More particulary, there is a race condition which may cause this
  //function to erroneously return true in that case.
  bool AllDone() const;
  

 private:
  class FeedThread : public IOneThread
  {
  public:
	  FeedThread(ThreadHandler *owner) : owner(owner) {}
  protected:
	  ThreadHandler *owner;

  	virtual bool OnDie() {
  		return true;
  	}

  	virtual bool OnDo(const string & msg) {
  		while(true) {
  			for(int i = 0; i < owner->GetNumThreads(); i++) {
  				if(owner->m_threads[i]->Done()) {
  					owner->m_threads[i]->Do(msg);
  					return true;
  				}
  			}
  			usleep(10000);
  		}
  	}

  	virtual bool OnInitialize(const string & msg) {
  		return true;
  	}
  };
  friend class FeedThread;

  FeedThread *feedThread;
  vector<IOneThread *> m_threads;
};




#endif //THREADHANDLER_H


