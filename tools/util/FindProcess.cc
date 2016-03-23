#include "base/FileParser.h"


#include <dirent.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "util/FindProcess.h"

#include <errno.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <stddef.h>


/*
int
readlink (const char *path, char *buf, size_t bufsize)
{
  struct stat statbuf;

  if (stat (path, &statbuf) >= 0)
    errno = EINVAL;
  return -1;
  }*/

int getProcessCount(const string & p_processname) {
#ifdef __APPLE_CC__
  printf("You are running MAC OS, checking for running applications is NOT supported.\n");
  return 0;
#endif


  DIR *dir_p;
  struct dirent *dir_entry_p;
  char dir_name[40];										// ??? buffer overrun potential
  char target_name[252];									// ??? buffer overrun potential
  int target_result;
  char exe_link[252];
  int errorcount;
  int result = 0xFFFFFFFF;
  
  errorcount=0;
  result=0;
  dir_p = opendir("/proc/"); 																// Open /proc/ directory

  int count = 0;

  while(NULL != (dir_entry_p = readdir(dir_p))) {											// Reading /proc/ entries
    if (strspn(dir_entry_p->d_name, "0123456789") == strlen(dir_entry_p->d_name)) {		// Checking for numbered directories 
      strcpy(dir_name, "/proc/");
      strcat(dir_name, dir_entry_p->d_name);
      strcat(dir_name, "/"); 															// Obtaining the full-path eg: /proc/24657/ 
      exe_link[0] = 0;
      strcat(exe_link, dir_name);
      strcat(exe_link, "exe");													 	// Getting 
      
      //printf("Checking %s\n", dir_name);
      strcat(dir_name, "status");
      
      FILE * pStat = fopen(dir_name, "r");
      if (pStat == NULL)
	continue;
      fclose(pStat);

      
      FlatFileParser parser;
      
      try {
	parser.Open(dir_name);
      }
      catch(...) {
	continue;
      }
      while (parser.ParseLine()) {
	if (parser.AsString(0) == "Name:") {
	  if (parser.AsString(1) == p_processname)
	    count++;
	}
      }
      
    }
  }
  closedir(dir_p);
  
  return count;
}

