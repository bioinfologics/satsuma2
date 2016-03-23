#ifndef XCORRDYNPROG_H
#define XCORRDYNPROG_H


#include "analysis/SequenceMatch.h"
#include "analysis/DNAVector.h"
#include <math.h>

typedef vec<char> qualvector;
typedef vec<qualvector> vecqualvector;

class MultiProtein;

#define XCINFINITY 999999999999999.;

class MDItem
{
public:
  MDItem() {
    m_qIndex = -1;
    m_score = XCINFINITY;
    m_prevIndex = -1;
    m_matchPenalty = XCINFINITY;
  }

  MDItem(int qIndex, double match) {
    m_qIndex = qIndex;
    m_score = XCINFINITY;
    m_prevIndex = -1;
    m_matchPenalty = match;
  }

  void Merge(int prev, double score, double penalty) {
    if (score + penalty < m_score) {
      m_prevIndex = prev;
      m_score = score + penalty;
    }
  }

  int GetIndex() const {return m_qIndex;}
  int GetPrevIndex() const {return m_prevIndex;}
  double GetScore() const {return m_score + m_matchPenalty;}

private:
  int m_qIndex;
  double m_score;
  double m_matchPenalty;
  int m_prevIndex;
};


class MDItemList
{
public:
  MDItemList() {

  }

  void Add(const MDItem & m) {
    m_nodes.push_back(m);
  }

  int GetCount() const {return m_nodes.isize();}
  MDItem & Get(int i) {return m_nodes[i];}


private:
  svec<MDItem> m_nodes;
};




class MatchDynProg
{
public:
  MatchDynProg(int targetStart, int targetLen);
  
  
  void AddMatch(const DNAVector & target, const DNAVector & query, const SingleMatch & m);
  int Merge(svec<int> & index);
  double PrettyPrint(const svec<int> & index, int start, const DNAVector & target, const DNAVector & query, int len);
  
private:
  svec<MDItemList> m_data;
  int m_targetStart;
};
 


class SigItem
{
 public:
  SigItem() {
    m_pos = 0;
    m_val = 0.;
  }

  void Set(int pos, double val) {
    m_pos = pos;
    m_val = val;
  }

  bool operator < (const SigItem & s) const {
    return m_val < s.m_val;
  }

  double Val() const {return m_val;}
  int Pos() const {return m_pos;}

  SigItem & operator = (const SigItem & s) {
    m_val = s.m_val;
    m_pos = s.m_pos;
    return *this;
  }

 private:
  double m_val;
  int m_pos;
};

class SignalFilter
{
 public:
  SignalFilter() {}

  int GetTopCount() const {return 15;}

  void Do(const svec<float> & signal) {
    int i;
    if (signal.isize() != m_items.isize())
      m_items.resize(signal.isize());
    for (int i=0; i<signal.isize(); i++) 
      m_items[i].Set(i-signal.isize()/2, signal[i]);

    Sort(m_items);

    /*
    int best = GetShift(0);
    int wiggle = 25;
    int k = m_items.isize()-2;
    for (i=m_items.isize()-2; i>=m_items.isize()-2-GetTopCount(); i--) {
      while (Bad(k, best, wiggle)) {     
	k--;
	if (k == 0)	  
	  break;
      }
     
      m_items[i] = m_items[k];
      k--;
      if (k <= 0)
	break;
	}*/
  }

  int GetShift(int i) const {
    return m_items[m_items.isize()-1-i].Pos();
  }
  


 private:
  bool Bad(int i, int pos, int max) const {
    int v = m_items[i].Pos();
    if (v - pos > max || pos - v > max)
      return true;
    return false;
  }

  svec<SigItem> m_items;
};


class XCDynProgLine
{
 public:
  XCDynProgLine() {
  }

  int Size() const {return m_score.isize();}

  void SetUp(const DNAVector & a, int shift, int size) {
    //cout << "Enter setup, size=" << size << endl;
    m_score.clear();
    m_letter.clear();
    m_back.clear();
    m_score.resize(size, 9999999999.);
    m_letter.resize(size, 0);
    m_back.resize(size, -1);
    //cout << "Done resize." << endl;

    bool first = true;
    for (int i=0; i<a.isize(); i++) {
      int x = i-shift;
      //cout << "i=" << i << " x=" << x << endl;
      if (x < 0 || x >= size) 
	continue;

      if (first) {
	m_score[x] = (double)x;
	first = false;
      }
      m_letter[x] = a[i];      
    }
    m_shift = shift;

    /*
    for (int j=0; j<m_letter.isize(); j++) {
      if (m_letter[j] == 0)
	cout << "-";
      else
	cout << m_letter[j];
    }
    cout << endl;*/

  }

  int Shift() const {return m_shift;}

  double Score(int i) const {return m_score[i];}
  char Letter(int i) const {return m_letter[i];}
  int Back(int i) const {return m_back[i];}

  void SetScore(int i, double s) {m_score[i] = s;}
  void SetBack(int i, int b) {m_back[i] = b;}


 private:
  int m_shift;
  svec<double> m_score;
  svec<char> m_letter;
  svec<int> m_back;
};



class XCDynProg
{
 public:
  XCDynProg() {
    m_thresh = 0.9;
    m_expect = 0.2;
  }

  void SetExpect(double d) {m_expect = d;}

  double Align(const MultiProtein & target, const MultiProtein &query, const SignalFilter & filter);

 private:

  bool FindBestBracket(int & lastFrame, int & m, int & firstFrame);

  double NormScore(double diff, int nodeDepth) {
    double var = sqrt((double)nodeDepth);
    double div = (double)nodeDepth-var;
    if (div <= 0.001)
      div = 0.001;
    return diff / div;
  }

  double m_thresh;
  double m_expect;

  svec<XCDynProgLine> m_matrix;


};



#endif
