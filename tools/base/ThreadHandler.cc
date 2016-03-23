#include "ThreadHandler.h"


void  * start_thread(void * str) {
  IOneThread * pMe = (IOneThread*)str;
  pMe->CallInit();
  return NULL;
}




bool IOneThread::Initialize(const string & msg, const string &name) {
  m_init = msg;
  this->name = name;
  //pthread_attr_t attr;
  //pthread_attr_init(&attr); /* initialize attr with default attributes */
  //pthread_attr_setscope(&attr, PTHREAD_SCOPE_SYSTEM);  /* system-wide contention */
  
  int iReturnValue1 = pthread_create(&m_iThreadId, NULL /*&attr*/, &start_thread, (void *)this);
  return true;
}


void IOneThread::CallInit() {
  OnInitialize(m_init);
  int i;
  
  do {
    bool bDo = false;
    string msg;   
    
    pthread_mutex_lock(&m_lock);
    if (m_msg.size() == 0) {
      bDo = false;
      m_bDo = false;
    } else {
      bDo = true;
      msg = m_msg[0];
      for (i=1; i<(int)m_msg.size(); i++)
	m_msg[i-1] = m_msg[i];
      m_msg.resize(m_msg.size()-1);
    }
    pthread_mutex_unlock(&m_lock);
    
    if (bDo) {
      OnDo(msg);
    } else {
      usleep(10000);
    }
    
  } while(!ShouldDie());
  OnDie();
  pthread_mutex_lock(&m_lock);
  m_bFinished = true;
  pthread_mutex_unlock(&m_lock);
}


//============================================================

ThreadHandler::ThreadHandler() {
	feedThread = NULL;
}

ThreadHandler::~ThreadHandler()
{
	if(feedThread != NULL) {
		while(!feedThread->Done())
			usleep(10000);
		feedThread->Die();
		while(!feedThread->Finished())
			usleep(10000);
		delete feedThread;
	}

	int i;

	//Wait until finished...
	int n = 0;
	do {
		n = 0;
		for (i=0; i<(int)m_threads.size(); i++) {
			if (!m_threads[i]->Done())
				n++;
		}
		if (n > 0)
			usleep(10000);

	} while (n > 0);

	for (i=0; i<(int)m_threads.size(); i++) {
		m_threads[i]->Die();
	}

	do {
		n = 0;
		for (i=0; i<(int)m_threads.size(); i++) {
			if (!m_threads[i]->Finished())
				n++;
		}
		if (n > 0)
			usleep(10000);
	} while(n >0);

	for (i=0; i<(int)m_threads.size(); i++) {
		delete m_threads[i];
	}
}


int ThreadHandler::AddThread(IOneThread * pThread, const string & init_msg, const string &name)
{
  m_threads.push_back(pThread);
  pThread->Initialize(init_msg, name);
  return (int)m_threads.size()-1;
}

bool ThreadHandler::AllDone() const
{
  if(feedThread != NULL && !feedThread->Done())
	  return false;
  for (int i=0; i<(int)m_threads.size(); i++) {
    if (!m_threads[i]->Done()) {
      return false;
    }
  }
  return true;
}

bool ThreadHandler::Feed(int threadID, const string & msg)
{
  return m_threads[threadID]->Do(msg);
}

void ThreadHandler::FeedAny(const string & msg) {
	if(feedThread == NULL) {
		feedThread = new FeedThread(this);
		feedThread->Initialize("", "Feed Thread");
	}
	feedThread->Do(msg);
}

 

