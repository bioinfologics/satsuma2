#ifndef RANDOMSTUFF_H_
#define RANDOMSTUFF_H_

#include <math.h>
#include "base/SVector.h"

inline long long RandomInt(long long n)
{    
  long long r =   ((long long)( rand( ) ) << 31 ) ^ (long long)( rand( ) );
  //cout << "r=" << r;
  r = r % n;
  //cout << " " << r << endl;
  return r;
}

inline double RandomFloat(double n)
{    
  long long r =   ((long long)( rand( ) ) << 31 ) ^ (long long)( rand( ) );
  //cout << "r=" << r;
  r = r % 0x7FFFFFFF;
  //cout << " " << r << endl;
  return n*((double)r)/(double)0x7FFFFFFF;
}




template<class T> 
void SRandomize(svec<T>& v)
{
  random_shuffle(v.begin(), v.end());
}

template<class T> 
void Randomize(vector<T>& v)
{
  random_shuffle(v.begin(), v.end());
}







#endif //RANDOMSTUFF_H_

