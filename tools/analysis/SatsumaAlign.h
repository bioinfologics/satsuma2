
#ifndef _SATSUMAALIGN_H_
#define _SATSUMAALIGN_H_

#include "base/SVector.h"

#include <string>

class SAlign
{
 public:
  SAlign() {
    m_startTarget = -1;
    m_endTarget = -1;
    m_startQuery = -1;
    m_endQuery = -1;
    m_direction = 0;
    m_ident = 0.;
  }


  void Set(const string & t, 
	   const string & q,
	   int startTarget,
	   int endTarget,
	   int startQuery,
	   int endQuery,
	   double ident,
	   const string & dir) {

    m_target = t;
    m_query = q;
    m_startTarget = startTarget;
    m_endTarget = endTarget;
    m_startQuery = startQuery;
    m_endQuery = endQuery;
    m_ident = ident;
    if (dir == "+" || dir == "1")
      m_direction = 1;
    if (dir == "-" || dir == "-1")
      m_direction = -1;

  }

  
  const string & Query() const {return m_query;}
  const string & Target() const {return m_target;}
  int StartTarget() const {return m_startTarget;}
  int EndTarget()   const {return m_endTarget;}
  int StartQuery()  const {return m_startQuery;}
  int EndQuery()    const {return m_endQuery;}
  int Direction()   const {return m_direction;}
  double Identity() const {return m_ident;}


 private:
  string m_query;
  string m_target;
  int m_startTarget;
  int m_endTarget;
  int m_startQuery;
  int m_endQuery;
  int m_direction;
  double m_ident;

};




class SAlignVec
{
 public:
  SAlignVec() {
    m_count = 0;
  }

  void Read(const string & fileName);
  void Write(const string & fileName);

  SAlign & operator[] (int i) {return m_aligns[i];}
  const SAlign & operator[] (int i) const {return m_aligns[i];}

  int isize() const {return m_count;}

  void clear() {
    m_count = 0;
    m_aligns.clear();
  }

 private:
  int m_count;
  svec<SAlign> m_aligns;
};


#endif //_SATSUMAALIGN_H_

