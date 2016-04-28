#include "SComm.h"

#include <stdio.h>

#include <stdlib.h>

#include <unistd.h>

#include <errno.h>
#include <string.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <arpa/inet.h>
#include <sys/wait.h>
#include <signal.h>
#include <netdb.h>

// max number of bytes we can get at once
#define MAXDATASIZE 1024


void w_error(const char * msg) {
  printf("ERROR: %s\n", msg);
  fflush(stdout);
}

class SSocketCommTransmitter : public SCommTransmitter
{
 public:
  SSocketCommTransmitter(int port = MYPORT);

  virtual ~SSocketCommTransmitter();

  virtual bool SendWait(const char * message);

private:
  /* listen on sock_fd, new connection on new_fd */
  int sockfd, new_fd;

  /* my address information */
  struct sockaddr_in my_addr;

  /* connector's address information */
  struct sockaddr_in their_addr;

  socklen_t sin_size;

  struct sigaction sa;

  int yes;
  int m_port;
};


class SSocketCommReceiver : public SCommReceiver
{
 public:
  SSocketCommReceiver(const char * serverName, int port = MYPORT);
  virtual ~SSocketCommReceiver() {}

  virtual bool Get(char * message, int bufSize);
  virtual bool IsFatal() {return m_bFatal;}

private:
  char m_serverName[512];
  int m_port;
  bool m_good;
  bool m_bFatal;
};


SCommTransmitter * GetTransmitter(int port)
{
  return new SSocketCommTransmitter(port);
}

SCommReceiver * GetReceiver(const char * serverName, int port)
{
  return new SSocketCommReceiver(serverName, port);
}

bool GetHostName(char * hostName, int bufSize)
{
  gethostname(hostName, bufSize);
  //printf("Local host: %s\n", hostName);
  return true;
}


/* how many pending connections queue will hold */
#define BACKLOG 10
void sigchld_handler(int s)
{  
  while(wait(NULL) > 0); 
}


//=========================================================================
SSocketCommTransmitter::SSocketCommTransmitter(int port)
{
  m_port = port;
  
  yes = 1;


  if ((sockfd = socket(AF_INET, SOCK_STREAM, 0)) == -1)
  {
    w_error("Server-socket() error lol!");    
    //exit(1);
  } else {
    //printf("Server-socket() sockfd is OK...\n");
  }

  if (setsockopt(sockfd, SOL_SOCKET, SO_REUSEADDR, &yes, sizeof(int)) == -1)
  {
    w_error("Server-setsockopt() error lol!");
    //exit(1);
  } else {
    //printf("Server-setsockopt is OK...\n");
  }
  
  
  /* host byte order */
  my_addr.sin_family = AF_INET;

  /* short, network byte order */
  my_addr.sin_port = htons(m_port);
  
  /* automatically fill with my IP */
  my_addr.sin_addr.s_addr = INADDR_ANY;

  //printf("Server-Using %s and port %d...\n", inet_ntoa(my_addr.sin_addr), m_port);

     

  /* zero the rest of the struct */
  memset(&(my_addr.sin_zero), '\0', 8);

  if(bind(sockfd, (struct sockaddr *)&my_addr, sizeof(struct sockaddr)) == -1)    
  {
    w_error("Server-bind() error");
    //exit(1);
  }
  else {
    //printf("Server-bind() is OK...\n");     
  }
  if(listen(sockfd, BACKLOG) == -1)
  {
    w_error("Server-listen() error");
    //exit(1);
  }
  //printf("Server-listen() is OK...Listening...\n");

  /* clean all the dead processes */
  sa.sa_handler = sigchld_handler;
  sigemptyset(&sa.sa_mask);
  sa.sa_flags = SA_RESTART;

  if(sigaction(SIGCHLD, &sa, NULL) == -1)
  {
    w_error("Server-sigaction() error");
    //exit(1);
  } else {
    //printf("Server-sigaction() is OK...\n");
  }
}

SSocketCommTransmitter::~SSocketCommTransmitter() 
{
  shutdown(sockfd, 2);
  close(sockfd);
}


bool SSocketCommTransmitter::SendWait(const char * message)
{
  int len = strlen(message) + 1;

  /* accept() loop */
  //while(1)
  //{
  sin_size = sizeof(struct sockaddr_in);
  
  if((new_fd = accept(sockfd, (struct sockaddr *)&their_addr, &sin_size)) == -1)
  {
    w_error("Server-accept() error");
    return false;      
  } else {
    //printf("Server-accept() is OK...\n");
  }
  //printf("Server-new socket, new_fd is OK...\n");
  
  //printf("Server: Got connection from %s\n", inet_ntoa(their_addr.sin_addr));
  
  /* this is the child process */
  if(!fork())
  {
    /* child doesn't need the listener */
    close(sockfd);

#ifdef __APPLE__
#define MSG_CONFIRM 0
#endif
    int ret = send(new_fd, message, len, MSG_CONFIRM);
    printf("Send returns: %d\n", ret);
    if (ret < 0) 
      w_error("Server-send() error lol!");
    close(new_fd);
    exit(0);
  } else {
    //printf("Server-send is OK...!\n");
  }
  /* parent doesn't need this*/
  
  close(new_fd);
  
  //printf("Server-new socket, new_fd closed successfully...\n");
    //}

  return true;
}
//=========================================================================



SSocketCommReceiver::SSocketCommReceiver(const char * serverName, int port)
{
  strcpy(m_serverName, serverName);
  m_port = port;
  m_good = true;
  m_bFatal = false;
}

bool SSocketCommReceiver::Get(char * message, int bufSize)
{
  
  int sockfd, numbytes;
  //char buf[MAXDATASIZE];
  struct hostent *he;
  
  // connector's address information
  
  struct sockaddr_in their_addr;
  
  
  // get the host info
  if((he=gethostbyname(m_serverName)) == NULL) {
    w_error("gethostbyname()");
    return false;
  } else {
    printf("Client-The remote host is: %s\n", m_serverName);
    fflush(stdout);
  }
  
  
  if((sockfd = socket(AF_INET, SOCK_STREAM, 0)) == -1)
  {
    w_error("socket()");
    return false;
  } else {
    printf("Client-The socket() sockfd is OK...\n");
    fflush(stdout);
  }
   
   
  // host byte order
  
  their_addr.sin_family = AF_INET;
  
  // short, network byte order
  printf("Server-Using %s and port %d...\n", m_serverName, m_port);
  
  their_addr.sin_port = htons(m_port);
  their_addr.sin_addr = *((struct in_addr *)he->h_addr);
  
  // zero the rest of the struct
  memset(&(their_addr.sin_zero), '\0', 8);     
  
  if(connect(sockfd, (struct sockaddr *)&their_addr, sizeof(struct sockaddr)) == -1)
  {
    w_error("connect()");
    close(sockfd);
    shutdown(sockfd, 2);
    return false;
  } else {
    printf("Client-The connect() is OK...\n");
    fflush(stdout);
  }
   

  struct timeval tv;
  //memset(&tv, sizeof(tv), 0);
  tv.tv_sec = 5;
  tv.tv_usec = 0;
  if (setsockopt(sockfd, SOL_SOCKET, SO_RCVTIMEO, (char *)&tv,  sizeof tv))
  {
    close(sockfd);
    w_error("setsockopt");
    return false;
  }


  if((numbytes = recv(sockfd, message, bufSize-1, 0)) == -1)
  {	
    w_error("recv()");
    m_bFatal = true;
    close(sockfd);
    return false;
  } else {
    //printf("numbytes=%d\n", numbytes);
    printf("Client-The recv() is OK...\n");
    fflush(stdout);
 }

  message[numbytes] = '\0';
  
  //printf("Client-Received: %s", message);     
  
  printf("Client-Closing sockfd\n");
  fflush(stdout);
  
  close(sockfd);
  shutdown(sockfd, 2);
  
  return true;
   
}
