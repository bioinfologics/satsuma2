#include "util/SysTime.h"
#include <unistd.h>
#include <strstream>
#include <iomanip>
#include <cstring>

using namespace std;

// this const was part of a class...
const int getTimeOfDayBufferLength = 27;

// this function was a static method in that class...
char* GetTimeOfDayInt(char* buffer, int bufferLength)
{
  if (bufferLength<getTimeOfDayBufferLength)
  {
    buffer[0] = '\0';
  }
  else
  {
    timeval tv;
    gettimeofday( &tv, NULL );
    strftime( buffer, 20,
      "%Y/%m/%d %H:%M:%S", localtime(&tv.tv_sec) );
    ostrstream ostr;
    ostr << ':' << setfill('0') << setw(6) << tv.tv_usec
      << ends;
    char* sp = ostr.str();
    strcat(buffer, sp);
    delete [] sp;
  }
  return buffer;
}

void GetTime(string & s) 
{
  char buffer[256];
  GetTimeOfDayInt(buffer, sizeof(buffer));
  s = buffer;
}

const string GetTimeStatic(bool bQuiet) 
{
  if (bQuiet) 
    return "";

  char buffer[256];
  GetTimeOfDayInt(buffer, sizeof(buffer));
  string s;
  s = buffer;
  return s;
}
