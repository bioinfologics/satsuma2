#ifndef _DPALIGN_H_
#define _DPALIGN_H_


#include "analysis/DNAVector.h"

class FullAlignment
{
 public:
  FullAlignment() 
    : m_offsetA(0), m_offsetB(0) {}

  private:
  string m_nameA;
  string m_nameB;

  svec<int> m_indexA;
  svec<int> m_indexB;
  svec<int> m_direction;
  svec<int> m_question;

  string m_strand;
 
  int m_offsetA, m_offsetB;

  void FlipOne(svec<int> & v, int qLen) {   
    for (int i=0; i<v.isize(); i++) {
      if (v[i] >= 0)
	v[i] = qLen - 1 - v[i];     
    }
  }


  public:

  const string & Strand() const {return m_strand;} 
  void SetStrand(const string & s) {m_strand = s;} 

  int Size() const {return m_indexA.isize();}

  int IndexA(int i) const {return (m_indexA[i]>-1 ? m_indexA[i]+m_offsetA : -1);}
  int IndexB(int i) const {return (m_indexB[i]>-1 ? m_indexB[i]+m_offsetB : -1);}
  int Direction(int i) const {return m_direction[i];}
  int Quest(int i) const {return m_question[i];}

  int OffsetA() const {return m_offsetA; }
  int OffsetB() const {return m_offsetB; }

  void resize(int n) {
    m_indexA.resize(n, -1);
    m_indexB.resize(n, -1);
    m_direction.resize(n, 0);
    m_question.resize(n, 0);
  } 


  void SetA(int i, int index) {m_indexA[i] = index;}
  void SetB(int i, int index) {m_indexB[i] = index;}
  void SetDirect(int i, int dir) {m_direction[i] = dir;}
  void SetQueston(int i, int q) {m_question[i] = q;}

  void SetOffsetA(int i) { m_offsetA = i; }
  void SetOffsetB(int i) { m_offsetB = i; }

  void Set(int i, int indexA, int indexB, int dir) 
  {
    m_indexA[i] = indexA;
    m_indexB[i] = indexB;
    m_direction[i] = dir;
  }

  void Push(int indexA, int indexB, int dir, int quest) {
    m_indexA.push_back(indexA);
    m_indexB.push_back(indexB);
    m_direction.push_back(dir);
    m_question.push_back(quest);
  }

  void SetNameA(const string &a) { m_nameA = a; }
  void SetNameB(const string &b) { m_nameB = b; }

  string GetNameA() const { return m_nameA; }
  string GetNameB() const { return m_nameB; }


  friend bool operator< (const FullAlignment &lhs, const FullAlignment &rhs)
  {
    if (lhs.m_nameA < rhs.m_nameA)
      return true;
    
    if ( lhs.m_nameA > rhs.m_nameA)
      return false;

    if ( lhs.IndexA(0) < lhs.IndexA(0) )
      return true;

    return false;

  }

  void WriteBinary(ofstream &oStrm);
  void ReadBinary(ifstream &iStrm);
  double Identity(const DNAVector & a, 
		  const DNAVector & b,
		  bool rcExplicit = false,
		  int offA = 0,
		  int offB = 0,
		  int printOff = 0);

  int PrettyPrint(const DNAVector & a, 
		  const DNAVector & b,
		  bool rcExplicit = false,
		  int offA = 0,
		  int offB = 0,
		  int printOff = 0);


  void Flip(int qLen) {
    FlipOne(m_indexB, qLen);
    for (int i=0; i<m_direction.isize(); i++)
      m_direction[i] = -m_direction[i];
  }



  void clear() {
    m_indexA.clear();
    m_indexB.clear();
    m_direction.clear();
    m_question.clear();

    m_nameA.erase();
    m_nameB.erase();
  }

};


struct orderFullAlignmentByNameB 
  : public binary_function< const FullAlignment, const FullAlignment, bool>
{
bool operator() (const FullAlignment &lhs, const FullAlignment &rhs) const
{
  return lhs.GetNameB() < rhs.GetNameB(); 
}
};

struct orderFullAlignmentByNameBandIndexB 
  : public binary_function< const FullAlignment, const FullAlignment, bool>
{
bool operator() (const FullAlignment &lhs, const FullAlignment &rhs) const
{
  if ( lhs.GetNameB() < rhs.GetNameB() )
    return true;

  if ( lhs.GetNameB() > rhs.GetNameB() )
    return false;

  if ( lhs.IndexB(0) < rhs.IndexB(0) )
    return true;

  return false;
}
};

struct orderFullAlignmentByNameBandNameA 
  : public binary_function< const FullAlignment, const FullAlignment, bool>
{
bool operator() (const FullAlignment &lhs, const FullAlignment &rhs) const
{
  if ( lhs.GetNameB() < rhs.GetNameB() )
    return true;

  if ( lhs.GetNameB() > rhs.GetNameB() )
    return false;

  if ( lhs.GetNameA() < rhs.GetNameA() )
    return true;

  return false;
}
};


struct orderFullAlignmentByNameBandNameAandIndexB
  : public binary_function< const FullAlignment, const FullAlignment, bool>
{
bool operator() (const FullAlignment &lhs, const FullAlignment &rhs) const
{
  if ( lhs.GetNameB() < rhs.GetNameB() )
    return true;

  if ( lhs.GetNameB() > rhs.GetNameB() )
    return false;

  if ( lhs.GetNameA() < rhs.GetNameA() )
    return true;

  if ( lhs.IndexB(0) < rhs.IndexB(0) )
    return true;


  return false;
}
};


class AlignsPrinter
{
 public:
  AlignsPrinter() {
    k = 0;
  }

  void PrettyPrint(const FullAlignment & align,
		   const DNAVector & a, 
		   const DNAVector & b,
		   int offA = 0,
		   int offB = 0,
		   bool bNew = true);

  
 private:
  int k;
  string seqA;
  string seqB;
  string match;
  string dir;

};



//==================================================================

#define FLIP_PENALTY 5.


/*
// The fw/rc version - but it doesn't work.
class DPAlignNode
{
 public:
  DPAlignNode() {
    m_backColFW = m_backColRC = -1;
    m_backRowFW = m_backRowRC -1;
    m_scoreFW = m_scoreRC = 9999999999.;
    //m_backRC = 0;
  }

  bool MinimizeScoreFW(double score, int backCol = -1, int backRow = -1) {
    if (score < m_scoreFW) {
      m_scoreFW = score;
      m_backColFW = backCol;
      m_backRowFW = backRow;
      //m_backRC = backRC;
      return true;
    }
    return false;
  }

  bool MinimizeScoreRC(double score, int backCol = -1, int backRow = -1) {
    if (score < m_scoreRC) {
      m_scoreRC = score;
      m_backColRC = backCol;
      m_backRowRC = backRow;
      //m_backRC = backRC;
      return true;
    }
    return false;
  }

  int GetBackColFW() const {return m_backColFW;}
  int GetBackRowFW() const {return m_backRowFW;}
  double GetScoreFW() const {return m_scoreFW;}
  int GetBackColRC() const {return m_backColRC;}
  int GetBackRowRC() const {return m_backRowRC;}
  double GetScoreRC() const {return m_scoreRC;}

  int GetBackCol(int dir) const {
    if (dir == 1)
      return m_backColFW;
    else
      return m_backColRC;
  }

  int GetBackRow(int dir) const {
   if (dir == 1)
     return m_backRowFW;
   else 
      return m_backRowRC;
  }


  double GetScore(int dir) const {
    if (dir == 1)
      return m_scoreFW;
    else
      return m_scoreRC;
  }

  double GetBestScore() const {
    if (m_scoreFW <= m_scoreRC)
      return m_scoreFW;
    else
      return m_scoreRC;
  }

  int IsBestRC() const {
    if (m_scoreFW <= m_scoreRC)
      return 1;
    else
      return -1;
  }


 private:
  int m_backColFW;
  int m_backRowFW;
  float m_scoreFW;

  int m_backColRC;
  int m_backRowRC;
  float m_scoreRC;
  //int m_backRC;
}; */



class DPAlignNode
{
 public:
  DPAlignNode() {
    m_backCol = -1;
    m_backRow = -1;
    m_score = 9999999999.;
  }

  bool MinimizeScore(double score, int backCol = -1, int backRow = -1) {
    if (score < m_score) {
      m_score = score;
      m_backCol = backCol;
      m_backRow = backRow;
      return true;
    }
    return false;
  }


  int GetBackCol() const {return m_backCol;}
  int GetBackRow() const {return m_backRow;}
  double GetScore() const {return m_score + m_reward;}

  int IsBestRC() const {return 1;}
 
  void SetReward(double d) {m_reward = d;}

 private:
  int m_backCol;
  int m_backRow;
  double m_score;
  double m_reward;
  
  //int m_backRC;
};

class DPAlignNodeColumn
{
 public:
  DPAlignNodeColumn(int y_off = 0) {
    m_y = y_off;
  }

  void Set(int y_off) {
    m_y = y_off;
  }

  void resize(int n) {
    m_col.resize(n);
  }

  int isize() const {
    return m_col.isize();
  }

  const DPAlignNode & operator [] (int i) const {return m_col[i];}
  DPAlignNode & operator [] (int i) {return m_col[i];}


  void SetScore(double d) {
    for (int i=0; i<m_col.isize(); i++) {
      m_col[i].MinimizeScore(d);
    }
  }
  void SetScore(double d, int i) {
    m_col[i].MinimizeScore(d);
  }

  bool Minimize(int y, double score, int backCol, int backRow) {
    int i = y - m_y;
    if (i < 0 || i >= m_col.isize())
      return false;
    
    if (y < backRow) {
      cout << "FATAL!!" << " y=" << y << " back=" << backRow << endl;
      //assert(0);
    }

    return m_col[i].MinimizeScore(score, backCol, backRow);
  }


  int GetY() const {return m_y;}

 private:
  svec<DPAlignNode> m_col;
  int m_y;
};


class TraceBackItem
{
 public:
  TraceBackItem(int row, int col, int dir, double score, int q) {
    m_score = score;
    m_row = row;
    m_col = col;
    m_dir = dir;
    m_quest = q;
  }
  
  
  int Row() const {return m_row;}
  int Col() const {return m_col;}
  int Dir() const {return m_dir;}
  int Quest() const {return m_quest;}
  
  
 private:
  double m_score;
  int m_row;
  int m_col;
  int m_dir;
  int m_quest;
};


//=====================================================================
class DPAligner
{
 public:
  DPAligner() {
    m_x_off = 0;
    m_rcScore = 35.;
    m_maxScore = 500.;
    m_maxSkip = (int)(m_maxScore + 0.5);
    m_maxLongSkip = 2000;
    m_minSkip = 50;
    m_band = 50.;
    m_bUseRewards = true;
  }

  void SetMaxSkip(int i) {
    m_maxLongSkip = i;
  }
  void UseRewardFunc(bool b) {
    m_bUseRewards = b;
  }

  void SetUpColumn(int x_off, int len);
  void SetUpRow(int x, int y_off, int len);

  void SetRect(int x_len, int y_len);

  void Align(FullAlignment & align,
	     const DNAVector & col, 
	     const DNAVector & row);


 private:

  void SetReward(const DNAVector & col, 
		 const DNAVector & row);

  int GetSkipLen(double scoreDiff) {
    double d = 1. - scoreDiff / m_band;
    if (d < 0)
      d = 0;
    int l = (int)(m_minSkip + (m_maxLongSkip - m_minSkip) * d);
    return l;
  }

  double RCPen() {
    return m_rcScore;
  }

  double RCPen(int rc, int rc2) {
    if (rc == 0 || rc2 == 0)
      return 0.;
    if (rc != rc2)
      return m_rcScore;
    return 0.;
  } 
  double SkipPen(int i1, int i2) {
    int diff = i2 - i1;
    if (diff < 0)
      diff = -diff;
    diff--;
    double s = 0.5 + diff;
    //double localMax = 25;
    double localMax = m_maxScore;
    if (s > localMax) {
      double diff = s - localMax;
      s = localMax + diff  * 0.01;
    }
    //cout << "i1=" << i1 << " i2=" << i2 << " pen=" << s << endl;
    return s;
  } 
  double Score(int x1, int x2, int y1, int y2, int rc) {
    int diffX = x2 - x1;
    int diffY = y2 - y1;

    if (diffX == 1) {
      if (rc == -1) {
	if (diffY == -1)
	  return 0.;
	else
	  return m_rcScore;
      } else {
	if (rc == 0)
	  return 0.;

	if (diffY == 1)
	  return 0.;
	else
	  return m_rcScore;
      }
    }
    if (diffX != 0 && diffY != 0)
      cout << "ERROR!!!" << endl;
    
    double s = 2. + (double)(diffX + diffY);
    if (s > m_maxScore)
      s = m_maxScore;
    return s;
  }

  void TraceBack(FullAlignment & align);

  double TraceBackScore(int row, int col,
			const DNAVector & colBases, 
			const DNAVector & rowBases);

  svec<DPAlignNodeColumn> m_columns;
  
  double m_band;
  int m_x_off;
  double m_rcScore;
  double m_maxScore;
  int m_maxSkip;
  int m_maxLongSkip;
  int m_minSkip;
  bool m_bUseRewards;
};


//==================================================================



#endif //_DPALIGN_H_



