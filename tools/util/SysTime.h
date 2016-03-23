#ifndef _SYSTIME_H_
#define _SYSTIME_H_

#include <string>
#include <iostream>

#include <sys/time.h>

using namespace std;

void GetTime(string & s);
const string GetTimeStatic(bool bQuiet = false);

#endif 
