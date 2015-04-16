#ifndef GRIDSEARCH_H_
#define GRIDSEARCH_H_


#include "analysis/SeqChunk.h"


class GridSequence
{
 public:
  GridSequence() {
    m_first = 0;
    m_last = 0;
  }

  void SetName(const string & name, int size) {
    m_name = name;
    m_size = size;
  }
  void SetBlocks(int first, int last) {
    m_first = first;
    m_last = last;
  }

  int First() const {return m_first;}
  int Last() const {return m_last;}
  const string & Name() {return m_name;}
  int Size() const {return m_size;}

 private:

  int m_first;
  int m_last;
  string m_name;
  int m_size;
};



enum SEARCH_STATE {
  SEARCH_UNKNOWN,
  SEARCH_SUBMITTED,
  SEARCH_DONE,
  SEARCH_MISS,
  SEARCH_HIT
};

class SearchCoordinates
{
 public:
  SearchCoordinates() {
    m_targetFirst = 0;
    m_targetLast = 0;
    m_queryFirst = 0;
    m_queryLast = 0;
  }

  SearchCoordinates(int tFirst, int tLast, int qFirst, int qLast) {
    m_targetFirst = tFirst;
    m_targetLast = tLast;
    m_queryFirst = qFirst;
    m_queryLast = qLast;
  }

  void Set(int tFirst, int tLast, int qFirst, int qLast) {
    m_targetFirst = tFirst;
    m_targetLast = tLast;
    m_queryFirst = qFirst;
    m_queryLast = qLast;
  }


  int TargetFirst() const {return m_targetFirst;}
  int TargetLast()  const {return m_targetLast;}
  int QueryFirst()  const {return m_queryFirst;}
  int QueryLast()   const {return m_queryLast;}


 private:
  int m_targetFirst;
  int m_targetLast;
  int m_queryFirst;
  int m_queryLast;
  

};


class SearchMatrix
{
 public:
  SearchMatrix() {
    m_x = -1;
    m_y = -1;
    m_blockX = -1;
    m_blockY = -1;

  }

  void Setup(int nTarget, int nQuery, int blockX, int blockY) {
    m_blockX = blockX;
    m_blockY = blockY;

    m_x = nTarget / blockX + 1;
    m_y = nQuery / blockY + 1;

    m_matrix.resize(m_x * m_y, SEARCH_UNKNOWN);
    m_coords.resize(m_x * m_y);
    m_counts.resize(m_x * m_y, 0);

    cout << "size x=" << m_x << " y=" << m_y << " size=" << m_coords.isize() << endl;

    int i, j;
    for (i=0; i<m_x; i++) {
      int lastX = (i + 1) * blockX - 1;
      if (lastX >= nTarget)
	lastX = nTarget - 1;
      for (j=0; j<m_y; j++) {
	int lastY =  (j + 1) * blockY - 1;
	if (lastY >= nQuery)
	  lastY = nQuery - 1;
	//cout << "i=" << i << " j=" << j << " index=" << Index(i, j) << endl;
	m_coords[Index(i, j)].Set(i * blockX, lastX, j * blockY, lastY); 
      }
    }

  }

  void ClearCounts() {
    for (int i=0; i<m_counts.isize(); i++) 
      m_counts[i] = 0;	
  }


  


  const SearchCoordinates & Coordinates(int x, int y) const {return m_coords[Index(x, y)];}

  SEARCH_STATE Get(int x, int y) const {return m_matrix[Index(x, y)];}
  int GetCount(int x, int y) const {return m_counts[Index(x, y)];}
  void AddCount(int x, int y, int c) {
    int i = Index(x, y);
    if (c > m_counts[i])
      m_counts[i] = c;
  }

  void Set(int x, int y, SEARCH_STATE s) {m_matrix[Index(x, y)] = s;}
  const char * Print(int x, int y) const {
    switch(m_matrix[Index(x, y)]) {
    case SEARCH_UNKNOWN:
      return "UNKNOWN";
    case SEARCH_SUBMITTED:
      return "SUBMITTED";
    case SEARCH_DONE:
      return "DONE";
    case SEARCH_MISS:
      return "MISS";
    case SEARCH_HIT:
      return "HIT";
    }
    return NULL;
  }

  void BlockByAbsolute(int & x, int & y, int targetBlock, int queryBlock) {
    x = targetBlock / m_blockX;
    y = queryBlock / m_blockY;
  }

  int TargetBlocks() const {return m_x;}
  int QueryBlocks() const {return m_y;}


 private:
  int Index(int x, int y) const {
    return y * m_x + x;
  }

  int m_x;
  int m_y;
  int m_blockX;
  int m_blockY;
  svec<SEARCH_STATE> m_matrix;
  svec<int> m_counts;

  svec<SearchCoordinates> m_coords;
  //vec<int> m_candidateCount;
  

};


class GridTarget
{
 public:
  GridTarget() {
    m_targetFirst = 0;
    m_targetLast = 0;
    m_queryFirst = 0;
    m_targetLast = 0;
    m_x = m_y = m_y2 = -1;
    m_count = 0;
    m_bFast = false;
  }

  GridTarget(int targetFirst, int targetLast, int queryFirst, int queryLast, int x, int y, int count) {
    Set(targetFirst, targetLast, queryFirst, queryLast, x, y);
    m_count = count;
    m_y2 = -1;
    m_bFast = false;
      
  }
  GridTarget(int targetFirst, int targetLast, int queryFirst, int queryLast, int x, int y, int y2, int count, bool bFast) {
    Set(targetFirst, targetLast, queryFirst, queryLast, x, y);
    m_count = count;
    m_y2 = y2;
    m_bFast = bFast;
  }

  void Set(int targetFirst, int targetLast, int queryFirst, int queryLast, int x, int y) {
    m_targetFirst = targetFirst;
    m_targetLast = targetLast;
    m_queryFirst = queryFirst;
    m_queryLast = queryLast;
    m_x = x;
    m_y = y;
  }

  bool IsFast() const {return m_bFast;}
  int TargetFirst() const {return m_targetFirst;}
  int TargetLast() const  {return m_targetLast;}
  int QueryFirst() const {return m_queryFirst;}
  int QueryLast() const {return m_queryLast;}
  int X() const {return m_x;}
  int Y() const {return m_y;}
  int Y2() const {return m_y2;}

  int GetCount() const {return m_count;}
  void SetCount(int c) {m_count = c;}

  bool operator <= (const GridTarget & t) const {
    return !operator > (t);
    //if (m_count != t.m_count)
    //  return (-m_count <= -t.m_count); 
    //if (m_x != t.m_x)
    //  return (m_x <= t.m_x); 
    //return (m_y <= t.m_y);
  }

  bool operator > (const GridTarget & t) const {
    if (m_count != t.m_count)
      return (-m_count > -t.m_count); 
    if (m_x != t.m_x)
      return (m_x > t.m_x); 
    return (m_y > t.m_y);
  }

  bool operator < (const GridTarget & t) const {
    if (m_count != t.m_count)
      return (-m_count < -t.m_count); 
    if (m_x != t.m_x)
      return (m_x < t.m_x); 
    return (m_y < t.m_y);
  }





 private:
  int m_targetFirst;
  int m_targetLast;
  int m_queryFirst;
  int m_queryLast;
  int m_x;
  int m_y;
  int m_y2;

  int m_count;

  bool m_bFast;

};


class RepeatTrackerItem
{
 public:
  RepeatTrackerItem() {}

  void SetSize(int size) {
    m_count.resize(size, 0);
  }

  bool IsRepeat(int start, int end) {
    int c = (start + end) / 2;
    //cout << "Asking " << start << " - " << end << " c=" << c << " count=" << m_count[c] << endl;
    if (m_count[c] > 1)
      return true;
    return false;
  }

  void SetRepeat(int start, int end) {
    int i;
    if (end >= m_count.isize()) {
      //cout << "ERROR: start= " << start << " end=" << end << " size=" << m_count.isize() << endl;
      end = m_count.isize() - 1;
    }
    if (start < 0)
      start = 0;
    for (i=start; i<end; i++)
      m_count[i]++;
    //cout << "Set repeat start=" << start << " end=" << end << endl;
  }

 private:
  int m_gran;
  svec<int> m_count;
};


class RepeatTracker
{
 public:
  RepeatTracker() {}

  void Setup(int tSeqs, int qSeqs) {
    m_tReps.resize(tSeqs);
    m_qReps.resize(qSeqs);
  }

  void SetTargetSize(int i, int size) {
    m_tReps[i].SetSize(size);
  }
  void SetQuerySize(int i, int size) {
    m_qReps[i].SetSize(size);
  }

  bool IsRepeat(int targetID, int tStart, int tEnd, int queryID, int qStart, int qEnd) {
    if (m_tReps[targetID].IsRepeat(tStart, tEnd))
      return true;
    if (m_qReps[queryID].IsRepeat(qStart, qEnd))
      return true;
    return false;
  }

  void SetRepeat(int targetID, int tStart, int tEnd, int queryID, int qStart, int qEnd) {
    m_tReps[targetID].SetRepeat(tStart, tEnd);
    m_qReps[queryID].SetRepeat(qStart, qEnd);
  }


 private:
  svec<RepeatTrackerItem> m_tReps;
  svec<RepeatTrackerItem> m_qReps;
};


//========================================================================
class GridSearch
{
 public:
  GridSearch(int size, int blocks);

  void SetUp(const vecDNAVector & target, 
	     const vecDNAVector & query);

  
  int NTargetChunks() const {return m_targetChunks.isize();}
  int NQueryChunks() const {return m_queryChunks.isize();}

  void ConsiderTargets(int targetID,
		       int targetStart,
		       int targetEnd,
		       int queryID,
		       int queryStart,
		       int querytEnd,
		       double ident);

  void SetUsed(int targetID,
	       int targetStart,
	       int targetEnd,
	       int queryID,
	       int queryStart,
	       int querytEnd);

  int CollectTargets(svec<GridTarget> & targets, int n, int nSeeds = 0);

  void CollectSeeds(svec<GridTarget> & targets, int n);

  void ClearTargetWeights() {
    m_matrix.ClearCounts();
    int i;
    for (i=0; i<m_vertical.isize(); i++) 
      m_vertical[i] = 0;
    for (i=0; i<m_horizontal.isize(); i++) {
      if (m_horizontal[i] > 0)
	m_horizontal[i] = 1;
    }
  }


  bool GetSelect(SeqChunkSelect & vertical, svec<int> & horizontal, int howMany); 

  const SeqChunk & TargetChunk(int i) const {return m_targetChunks[i];}
  const SeqChunk & QueryChunk(int i)  const {return m_queryChunks[i];}

  

 private:

  svec<SeqChunk> m_targetChunks;
  svec<SeqChunk> m_queryChunks;
  int m_size;

  svec<GridSequence> m_targetSeq;
  svec<GridSequence> m_querySeq;
  RepeatTracker m_repTrack;

  int m_blocks;

  void CoordsToBlocks(int & targetBlock,
		      int & queryBlock,
		      int target,
		      int startTarget,
		      int query,
		      int startQuery);
  
  SearchMatrix m_matrix;

  svec<int> m_vertical;
  svec<int> m_horizontal;
  svec<int> m_allNs;

  svec<GridTarget> m_targets;

};


#endif //GRIDSEARCH_H_

