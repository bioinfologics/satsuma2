

#ifndef _SEQUENCEMATCH_H_
#define _SEQUENCEMATCH_H_

#include "../base/SVector.h"
#include <string>
#include <iostream>

// query is rc, not the target!
class CMReadFileStream;
class CMWriteFileStream;




class SingleMatch
{
 public:
  SingleMatch() {
    m_targetID = -1;
    m_queryID = -1;
    m_queryLen = 0;
    m_startTarget = 0;
    m_startQuery = 0;
    m_length = 0;
    m_bRC = false;
    m_matches = 0;
    m_prob = 0.;
    m_ident = 0.;
  }

  virtual ~SingleMatch() {}

  void Print() const {
    cout << "target=" << m_targetID << " @" << m_startTarget << " len=" << m_length << " query=";
    cout << m_queryID << " @" << m_startQuery << " prob=" << m_prob;
    if (m_bRC) {
      cout << " -";
    } else {
      cout << " +";
    }
    cout << endl;
  }


  void SetProbability(double p) {
    m_prob = p;
  }

  void SetQueryTargetID(int queryID, int targetID, int qLen) {
    m_queryID = queryID;
    m_targetID = targetID;
    m_queryLen = qLen;
  }

  
  void AddMatches(double d) {m_matches += d;}

  void SetPos(int queryStart, 
	      int targetStart, 
	      int len, 
	      bool rc) {
    m_startTarget = targetStart;
    m_length = len;   
    m_startQuery = queryStart;    
    m_bRC = rc;
  }

  void SetIdentity(double id) {
    m_ident = id;
  }

  int GetTargetID() const {return m_targetID;}
  int GetQueryID()  const {return m_queryID;}
 
  int GetStartTarget() const {return m_startTarget;}
  int GetStartQuery()  const {return m_startQuery;}
  int GetLength()      const {return m_length;}
  bool IsRC()          const {return m_bRC;}
  double GetMatches()  const {return m_matches;}
  double GetProbability()  const {return m_prob;}
  double GetIdentity()     const {return m_ident;}

  void SetLength(int l) {m_length = l;}

  void ReadAppend(CMReadFileStream & f, int version = 0, int app = 0);
  void Read(CMReadFileStream & f, int version = 0);
  void Write(CMWriteFileStream & f);

  virtual bool operator < (const SingleMatch &s) const {
    if (m_targetID != s.m_targetID)
      return (m_targetID < s.m_targetID);
    if (m_queryID != s.m_queryID)
      return (m_queryID < s.m_queryID);
    if (m_bRC != s.m_bRC) {
      if (m_bRC == false)
	return true;
      else
	return false;
    }

    if (m_startTarget == s.m_startTarget)
      return (m_length < s.m_length);

    return (m_startTarget < s.m_startTarget);
  }

  bool operator == (const SingleMatch &s) const {
    if (m_targetID != s.m_targetID)
      return false;
    if (m_queryID != s.m_queryID)
      return false;
    if (m_bRC != s.m_bRC) 
      return false;
    
    if (m_startTarget != s.m_startTarget)
      return false;

    if (m_startQuery != s.m_startQuery)
      return false;

    return true;
  }

  bool Close(const SingleMatch &s, int size) const {
    if (m_targetID != s.m_targetID)
      return false;
    if (m_queryID != s.m_queryID)
      return false;
    if (m_bRC != s.m_bRC) 
      return false;
    
    int thresh = size;
    int d = m_startTarget - s.m_startTarget;
    if (d < -thresh || d > thresh)
      return false;

    d = m_startQuery - s.m_startQuery;
    if (d < -thresh || d > thresh)
      return false;

    return true;
  }

 private:
  int m_targetID;
  int m_queryID;
  int m_queryLen;
  int m_startTarget;
  int m_startQuery;
  int m_length;
  bool m_bRC;
  double m_matches;
  double m_prob;
  double m_ident;
};




class MultiMatches
{
 public:
  MultiMatches() {
    m_count = 0;
    m_currVer = 3;
  }

  void SetNames(const svec<string> & qNames, const svec<string> & tNames);

  const string & GetQueryName(int i) const {return m_queryNames[i];}
  const string & GetTargetName(int i) const {return m_targetNames[i];}

  void SetTargetSize(int i, int size);
  void SetQuerySize(int i, int size);
  void SetTargetName(int i, const string & n);
  void SetQueryName(int i, const string & n);

  int GetTargetSize(int i) const {return m_targetSize[i];}
  int GetQuerySize(int i)  const {return m_querySize[i];}

  int GetTargetCount() const {return m_targetSize.isize();}
  int GetQueryCount()  const {return m_querySize.isize();}

  void SetCounts(int target, int query);


  void AddMatch(const SingleMatch & m) {
    //cout << "Added match." << endl;
    if (m_count >= m_matches.isize()) {
      m_matches.resize(m_count + 65536 * 2);
    }
    m_matches[m_count] = m;
    m_count++;
  }

  void Read(const string & file);
  void ReadSummary(const string & file);
  void MergeRead(const string & file);
  void Write(const string file);
  
  void MergeReadAppend(const string & file);

  int GetMatchCount() const {return m_count;}
  const SingleMatch & GetMatch(int i) const {return m_matches[i];}
  SingleMatch & GetMatchDirect(int i) {return m_matches[i];}

  void Sort() {
    m_matches.resize(m_count);
    ::Sort(m_matches);
  }


  void ClearMatches() {
    m_matches.clear();
    m_count = 0;
  }

  void Merge(const MultiMatches & m);


  void Clear() {
    m_count = 0;
    m_targetNames.clear();
    m_queryNames.clear();
    m_matches.clear();

    m_targetSize.clear();
    m_querySize.clear();
  }

  void Collapse();

 private:
  int Update(svec<string> & names, svec<int> & size, const string & n, int s);

  svec<string> m_targetNames;
  svec<string> m_queryNames;
  svec<SingleMatch> m_matches;
  int m_count;

  svec<int> m_targetSize;
  svec<int> m_querySize;
  int m_currVer;

};




#endif //_SEQUENCEMATCH_H_

