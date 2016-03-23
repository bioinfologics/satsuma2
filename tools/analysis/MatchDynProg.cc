#ifndef FORCE_DEBUG
#define NDEBUG
#endif


#include <string>
#include <math.h>
#include "../base/CommandLineParser.h"
//#include "analysis/CrossCorr.h"
#include "SequenceMatch.h"
#include "MatchDynProg.h"


#define PRETTY_INFINITE 999999999999999.;

#define RC_PENALTY 20.
#define MAX_PENALTY 20


typedef svec<char> qualvector;
typedef svec<qualvector> vecqualvector;


double MatchPenalty(double prob)
{
  double p = - 200 * log(prob);
  //cout << "prob=" << prob << " pen=" << p << endl;
  return p;
}

double GetMaxPenalty() 
{
  return MAX_PENALTY;
}

double GetHighMaxPenalty() 
{
  return 150;
}

double TransPenalty(int startT1, int startQ1, bool rc1, int startT2, int startQ2, bool rc2) 
{
  double pen = 0.;

  double expect = (double)(startT2 - startT1);
  double observe = (double)(startQ2 - startQ1);

  double sigma = sqrt(expect)/10.;

  sigma -= 100;
  if (sigma < 0)
    sigma = 0;

  if (rc1 != rc2) { 
    pen = RC_PENALTY;
  }

  if (rc1 == true && rc2 == true) {
    observe = -observe;
  }

  double diff = observe - expect;
  if (diff < 0.)
    diff = -diff;

  double flat = 0.01;
  if (diff == 1)
    flat = 0;

  if (expect < 300.)
    expect = 300;
  double p = diff / expect + sigma + flat;

  if (p > MAX_PENALTY)
    p = MAX_PENALTY;

  /*
  cout << "expect=" << expect << " observe=" << observe;
  if (rc1)
    cout << " rc1";
  else 
    cout << " fw1";
  if (rc2)
    cout << " rc2";
  else 
    cout << " fw2";
  
  cout << " pen+p=" << pen+p << endl;
  */

  return pen + p;
}

double GetRepeatScore(int v) {
  //    v += 2;
  if (v <= 1)
    return 10.;
  if (v == 2)
    return 5.;
  if (v == 3)
    return 2.5;
  if (v == 4)
    return 1.2;
  if (v >= 25)
    return 0.000001;
  if (v >= 10)
    return 0.001;

  return 0.5;
}


class SingleMatchDP : public SingleMatch
{

public:
  SingleMatchDP() {
    m_score = PRETTY_INFINITE;
    m_pen = -1.;
    m_back = -1;
    m_count = 0;
    m_repScore = 100.;
  }

  virtual bool operator < (const SingleMatch &s) const {
    return (GetStartTarget() < s.GetStartTarget());
  }

  void Update(double score, int back) {
    if (score < m_score) {
      m_score = score;
      m_back = back;
    }
  }

  void SetRepeatScore(double d) {
    m_repScore = d;
  }

  double GetRepeatScore() const {return m_repScore;}

  int GetBack() const {return m_back;}

  double GetScore() {
    if (m_pen < 0.)
      m_pen = MatchPenalty(GetProbability());
    return m_score + m_pen;
  }

  SingleMatchDP & operator = (const SingleMatch & m) {
    ((SingleMatch*)this)->operator =(m);
    return *this;
  }

private:
  double m_score;
  double m_pen;
  int m_back;

  int m_count;
  double m_repScore;
};


class MatchDynProg
{
public:
  MatchDynProg(int size) {
    m_matches.resize(size);
    m_laDist = 250000;
    m_minKeepLen = 220;
    m_minKeepIdent = 0.57;
  }

  void SetKeep(int minLen, double minIdent) {
    m_minKeepLen = minLen;
    m_minKeepIdent = minIdent;
  }

  void Set(const SingleMatch & m, int i, double d) {
    m_matches[i] = m;
    m_matches[i].SetRepeatScore(d);
  }

  void Close(int k) {
    m_matches.resize(k);
  }

  void Chain(svec<SingleMatch> & out);

private:
  svec<SingleMatchDP> m_matches;
  int m_laDist;
  int m_minKeepLen;
  double m_minKeepIdent;

};


void MatchDynProg::Chain(svec<SingleMatch> & out)
{
  cout << "Sorting..." << endl;
  Sort(m_matches);
  cout << "# of matches: " << m_matches.isize() << endl;
  if (m_matches.isize() == 0)
    return;


  int i, j;
  int la_limit = 2000;
  m_matches[0].Update(0., -1);
  for (i=0; i<m_matches.isize(); i++) {
    SingleMatchDP & one = m_matches[i];

    if (i % 1000 == 0)
      cout << "Processed: " << i << endl;

    double skipScore = 0.;
    int fed = 0;
    for (j=i+1; j<m_matches.isize(); j++) {
      SingleMatchDP & two = m_matches[j];
      if (two.GetStartTarget() - one.GetStartTarget() > m_laDist
	  && fed > 10) {
	//cout << "Forward break: " << j - i << endl;
	break;
      }
      if (j-i > la_limit)
	break;

      double trans = GetHighMaxPenalty();
      fed++;
      if (one.GetQueryID() == two.GetQueryID()) {
	trans = TransPenalty(one.GetStartTarget(), 
			     one.GetStartQuery(), 
			     one.IsRC(), 
			     two.GetStartTarget(), 
			     two.GetStartQuery(), 
			     two.IsRC());
      }



      trans += skipScore;
   
      //cout << "Dist=" << j-i << " score=" << skipScore << endl;
      
      skipScore += two.GetRepeatScore();

      //if (one.GetStartTarget() < 200000 && two.GetQueryID() != 0)
      //continue;
      two.Update(one.GetScore() + trans, i);

    }
  }
  //cout << "Dyn prog done, traceback." << endl;

  i = m_matches.isize()-1;
  while (i >= 0) {
    out.push_back(m_matches[i]);
    //cout << "Pos=" << m_matches[i].GetStartTarget()  << " score=" << m_matches[i].GetScore() << " query=" << m_matches[i].GetQueryID();
    //if (m_matches[i].IsRC())
    //  cout << " rc" << endl;
    //else
    //  cout << " fw" << endl;

    i = m_matches[i].GetBack();    
    //=============================================
    //m_matches[i].SetQueryTargetID(-1, -1, 0);
  }

  //cout << "NOT Adding forced keepers." << endl;
  int keep = 0;

  //=============================================
  //for (i=0; i<m_matches.isize(); i++) {
  //// We got this one already!
  //  if (m_matches[i].GetTargetID() == -1)
  //  continue;
  // if (m_matches[i].GetLength() > m_minKeepLen) {
  //   out.push_back(m_matches[i]);
  //   keep++;
  //   continue;
  // }
  // if (m_matches[i].GetIdentity() > m_minKeepIdent) {
  //   out.push_back(m_matches[i]);
  //   keep++;
  // }
  //}
  //cout << "Forced in " << keep << " matches." << endl;
  Sort(out);
}

bool RunMatchDynProgMult(MultiMatches & out, const MultiMatches & in)
{
  MultiMatches tmp;
  MultiMatches rest;
  
  RunMatchDynProg(tmp, in);

  tmp.Sort();

  int i, j;

  rest.SetCounts(in.GetTargetCount(), in.GetQueryCount());
  for (i=0; i<in.GetTargetCount(); i++) {
    rest.SetTargetSize(i, in.GetTargetSize(i));
    rest.SetTargetName(i, in.GetTargetName(i));
   }
  for (i=0; i<in.GetQueryCount(); i++) {
    rest.SetQuerySize(i, in.GetQuerySize(i));
    rest.SetQueryName(i, in.GetQueryName(i));
  }
  
  i = 0;
  j = 0;

  cout << "Removing dups" << endl;
  int skipped = 0;
  
 
  for (j=0; j<in.GetTargetCount(); j++) {

    svec<int> qHits;
    svec<int> qHitsHi;
    qHits.resize(in.GetQueryCount(), 0);
    qHitsHi.resize(in.GetQueryCount(), 0);

    for (i=0; i<tmp.GetMatchCount(); i++) {
      if (tmp.GetMatch(i).GetTargetID() != j)
	continue;
      int q = tmp.GetMatch(i).GetQueryID();
      if (q >= 0) {
	qHits[q]++;
	int pos = tmp.GetMatch(i).GetStartQuery();
	if (tmp.GetMatch(i).GetLength() > 20 && pos > qHitsHi[q])
	  qHitsHi[q] = pos;
      }
    }

    int bestQuery = -1;
    int bestQueryMax = -1;
    int maxHits = 0;
    

    for (i=0; i<qHits.isize(); i++) {
      if (qHits[i] > maxHits) {
	maxHits = qHits[i];
	bestQuery = i;
	bestQueryMax = qHitsHi[i];
      }      
    }
    cout << "Target: " << j << " best query: " << bestQuery << " Hits: " << maxHits << endl;
    cout << "  hi: " << bestQueryMax << endl;
  
    
    int scale = 500000;
    int len = 0;
    //int lo = 0;
    //int hi = 0;
    svec<int> numHitsScale;
    
    if (bestQuery != -1) {
      //numHits.resize(
      len = in.GetQuerySize(bestQuery);
      numHitsScale.resize(len/scale+1, 0);
      
      for (i=0; i<tmp.GetMatchCount(); i++) {
	if (tmp.GetMatch(i).GetTargetID() != j)
	  continue;
	int q = tmp.GetMatch(i).GetQueryID();
	if (q != bestQuery)
	  continue;
	int pos = tmp.GetMatch(i).GetStartQuery();
	if (tmp.GetMatch(i).IsRC())
	  pos = len - pos;
	numHitsScale[pos/scale]++;
      }
    }
    for (i=0; i<in.GetMatchCount(); i++) {
      if (in.GetMatch(i).GetTargetID() != j)
	continue;
      if (in.GetMatch(i).GetQueryID() != bestQuery) {
	rest.AddMatch(in.GetMatch(i));
      } else {
	/*
	int pos = in.GetMatch(i).GetStartQuery();
	if (in.GetMatch(i).IsRC())
	  pos = len - pos;
	if (bestQuery != -1 && numHitsScale[pos/scale] < 4)
	  rest.AddMatch(in.GetMatch(i));
	*/
      }
    }
  }
    
  /* 
  for (i=0; i<tmp.GetMatchCount(); i++) {
    if (i % 10000 == 0)
      cout << i << endl;
    while (j<in.GetMatchCount()) {
      if (in.GetMatch(j) == tmp.GetMatch(i)) {
	skipped ++;

	break;
      } 

      if (in.GetMatch(j).Close(tmp.GetMatch(i), 15000000)) {
	skipped++;
      } else {	
	rest.AddMatch(in.GetMatch(j));
      }
      j++;
      
    }

    }*/

    //cout << "done, skipped: " << skipped << endl;

  /*
  for (i=0; i<in.GetMatchCount(); i++) {
    if (i % 10000 == 0)
      cout << i << endl;
    bool bBad = false;
    for (j=0; j<tmp.GetMatchCount(); j++) {
      if (in.GetMatch(j).Close(tmp.GetMatch(i))) {
	bBad = true;
	break;
      }    
    }
    if (!bBad)
      rest.AddMatch(in.GetMatch(i));
     
  }
  */

  //out = rest;
  //return true;

  RunMatchDynProg(out, rest);

  for (i=0; i<tmp.GetMatchCount(); i++) 
    out.AddMatch(tmp.GetMatch(i));

  out.Sort();
  return true;
}

bool RunMatchDynProg(MultiMatches & out, const MultiMatches & in, bool bSelf)
{

  int i, j, k;

  cout << "Chaining (inline)..." << endl;

  int n = in.GetTargetCount();

  vecqualvector mult_target;
  vecqualvector mult_query;


  mult_target.resize(in.GetTargetCount());
  mult_query.resize(in.GetQueryCount());

  out.SetCounts(in.GetTargetCount(), in.GetQueryCount());
  for (i=0; i<in.GetTargetCount(); i++) {
    out.SetTargetSize(i, in.GetTargetSize(i));
    out.SetTargetName(i, in.GetTargetName(i));
    qualvector & q = mult_target[i];
    q.resize(in.GetTargetSize(i));
    for (j=0; j<in.GetTargetSize(i); j++)
      q[j] = 0;
  }
  for (i=0; i<in.GetQueryCount(); i++) {
    out.SetQuerySize(i, in.GetQuerySize(i));
    out.SetQueryName(i, in.GetQueryName(i));
    qualvector & q = mult_query[i];
    q.resize(in.GetQuerySize(i));
    for (j=0; j<in.GetQuerySize(i); j++)
      q[j] = 0;
  }

  cout << "Filling out repeat lists..." << endl;


  SingleMatch last;
  
  for (i=0; i<in.GetMatchCount(); i++) {
    const SingleMatch & m = in.GetMatch(i);

    if (m.GetQueryID() < 0 /*|| m.GetStartTarget() > 17000000*/) {
      //cout << "ERROR!" << endl;
      //cout << "target=" << m.GetTargetID() << " query=" << m.GetQueryID() << " t=" << m.GetStartTarget() << " q=" << m.GetStartQuery() << " len=" << m.GetLength() << endl;
     continue;
    }
    

    int startT = m.GetStartTarget();
    int startQ = m.GetStartQuery();

    if (m.GetTargetID() == last.GetTargetID() && m.GetQueryID() == last.GetQueryID() &&
	m.GetStartTarget() <= last.GetStartTarget() + last.GetLength() &&
	m.GetStartQuery() <= last.GetStartQuery() + last.GetLength() &&
	m.GetStartTarget() > last.GetStartTarget() &&
	m.GetStartQuery() > last.GetStartQuery()) {
      //last = m;
      //continue;
      startT = last.GetStartTarget() + last.GetLength();
      startQ = last.GetStartQuery() + last.GetLength();
      //cout << "Found..." << endl;
    }

    last = m;
    
    qualvector & t = mult_target[m.GetTargetID()];
    qualvector & q = mult_query[m.GetQueryID()];
    for (j=startT; j<m.GetStartTarget() + m.GetLength(); j++) {
      if (j>0 && j<t.isize()) {
	if (t[j] < 100)
	  t[j]++;
      }
    }
    for (j=startQ; j<m.GetStartQuery() + m.GetLength(); j++) {
      if (j>0 && j<q.isize()) {
	if (q[j] < 100)
	  q[j]++;
      }
    }      
  }

  cout << "Dynprog'ing..." << endl;


  svec<int> firstI, lastI;
  firstI.resize(n, in.GetMatchCount()+1);
  lastI.resize(n, -1);

 
  for (i=0; i<in.GetMatchCount(); i++) {
    const SingleMatch & m = in.GetMatch(i);    
    
    int id = m.GetTargetID();
    if (i < firstI[id])
      firstI[id] = i;
    if (i > lastI[id])
      lastI[id] = i;
  }
  
  
  // Less stupid way of doing this
  for (j=0; j<n; j++) {
    //cout << "Target " << j << " (of " << n << ")" << endl;
    int first = firstI[j];
    int last = lastI[j];

    //cout << "First=" << firstI[j] << " Last=" << lastI[j] << endl;

    if (last == -1)
      first = 0;

    /*
    for (i=0; i<in.GetMatchCount(); i++) {
      const SingleMatch & m = in.GetMatch(i);
      //if (m.GetQueryID() < 0)
      //continue;
      if (m.GetTargetID() != j)
	continue;
      if (first == -1)
	first = i;
      last = i;
    }
    */


    MatchDynProg dp(last - first + 1);

    int k = 0;

    dp.SetKeep(2220, 0.57);

    qualvector & t = mult_target[j];
 
    if (first == -1) {
      first = 0;
      last = -1;
    }
      
    //cout << "Collecting candidates." << endl;
    for (i=first; i<=last; i++) {
      //cout << "i=" << i << endl;
      //cout << "queryID=" << in.GetMatch(i).GetQueryID() << " Target id=" << in.GetMatch(i).GetTargetID() << endl;
      qualvector & q = mult_query[in.GetMatch(i).GetQueryID()];
      const SingleMatch & m = in.GetMatch(i);
 
      if (bSelf && m.GetQueryID() == m.GetTargetID())
	continue;

      /*if (m.GetQueryID() ==2  && m.GetStartTarget() > 17000000) {
	cout << "ERROR!" << endl;
	cout << "target=" << m.GetTargetID() << " query=" << m.GetQueryID() << " t=" << m.GetStartTarget() << " q=" << m.GetStartQuery() << " len=" << m.GetLength() << endl;
	continue;
      }*/
 


      //cout << "Rep score." << endl;
      //cout << "Size t: " << t.isize() << " " << m.GetStartTarget() << " " << m.GetLength()/2 << endl;
      //cout << "Size q: " << q.isize() << " " << m.GetStartQuery() << " " << m.GetLength()/2 << endl;
      double repScore1 = GetRepeatScore(t[m.GetStartTarget() + m.GetLength()/2]);
      double repScore2 = GetRepeatScore(q[m.GetStartQuery() + m.GetLength()/2]);
      if (repScore2 < repScore1)
	repScore1 = repScore2;
      //cout << "done, it's " << repScore1 << endl;

      if (repScore1 > 0.0001) {
	dp.Set(m, k, repScore1);
	k++;
      } else {
	//cout << "Skipping!" << endl;      
      }
    }

    //cout << "Done, chaining..." << endl;
    dp.Close(k);

    svec<SingleMatch> result;
    dp.Chain(result);
    

    for (i=0; i<result.isize(); i++) {
      const SingleMatch & m = result[i];
      out.AddMatch(m);
      //cout << "query id: " << m.GetQueryID() << " start target: " << m.GetStartTarget() << " start query: " << m.GetStartQuery();
      //if (m.IsRC()) {
      //	cout << " rc";
      //} else {
      //	cout << " fw";
      //}
      //cout << endl;
    }
  }

  return true;
}

