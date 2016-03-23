#ifndef _SCOMM_H_
#define _SCOMM_H_

/* the port users will be connecting to */
#define MYPORT 3491

class SCommTransmitter
{
 public:
  virtual ~SCommTransmitter() {}

  virtual bool SendWait(const char * message) = 0;
  
};


class SCommReceiver
{
 public:
  virtual ~SCommReceiver() {}

  // Make sure there is enough space in the buffer!
  virtual bool Get(char * message, int bufSize) = 0;
  virtual bool IsFatal() = 0;
};

// This is how to get an actual instance
SCommTransmitter * GetTransmitter(int port = MYPORT);
SCommReceiver * GetReceiver(const char * serverName, int port = MYPORT);

bool GetHostName(char * host, int bufSize);

#endif 


