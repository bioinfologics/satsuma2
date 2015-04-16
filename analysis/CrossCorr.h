#ifndef _CROSSCORR_H_
#define _CROSSCORR_H_

//#include "base/SVector.h"
#include "extern/RealFFT/FFTReal.h"
#include "analysis/DNAVector.h"




class CCSignal
{
 public:
  CCSignal();
  virtual ~CCSignal() {}

  virtual void SetSequence(const DNAVector & b, int size);

  virtual int GetFullSize() const {return GetA().size();}
  
  virtual void SetSize(int n);

  virtual void Smooth();


  virtual int GetCount() const {return 4;}
  virtual const svec<float> & Get(int i) const {
    switch(i) {
    case 0:
      return m_A;
      break;
    case 1:
      return m_C;
      break;
    case 2:
      return m_G;
      break;
    case 3:
      return m_T;
      break;
    }
    return m_A;
  }



 protected:
  svec<float> & A() {return m_A;};
  svec<float> & C() {return m_C;};
  svec<float> & G() {return m_G;};
  svec<float> & T() {return m_T;};

  const svec<float> & GetA() const {return m_A;};
  const svec<float> & GetC() const {return m_C;};
  const svec<float> & GetG() const {return m_G;};
  const svec<float> & GetT() const {return m_T;};
  


  void SeqToPCM(svec<float> & out, const DNAVector & in, char nuke); 
  void ComputeEntropy(const DNAVector & in);

  svec<float> m_A;
  svec<float> m_C;
  svec<float> m_G;
  svec<float> m_T;

  svec<float> m_entropy;
  int m_points;


};


inline double CCScore(const CCSignal & a, const CCSignal & b, int apos, int bpos)
{
  int i;
  double sum = 0;
  for (i=0; i<a.GetCount(); i++) {
    const svec<float> & af = a.Get(i);
    const svec<float> & bf = b.Get(i);
  
    double one = af[apos];
    double two = bf[bpos];
    double dot = one*two;
    sum += dot;
    //cout << "blah, i=" << i << " a=" << one << " b=" << two << endl;

     
  }
  //cout << "sum=" << sum << endl;
  return sum;
}

typedef svec<float> vecFloat;

class CCSignalProtein : public CCSignal
{
 public:
  CCSignalProtein() : CCSignal() {}
  virtual ~CCSignalProtein() {}

  virtual void SetSequence(const DNAVector & b, int size);

  virtual int GetFullSize() const {return m_aa[0].size();}
  
  virtual void SetSize(int n);

  virtual int GetCount() const {return 21;}
  virtual const svec<float> & Get(int i) const {
    return m_aa[i];
  }

 private:
  vecFloat m_aa[21];
};


class CCSignalWithCodons : public CCSignal
{
 public:
  CCSignalWithCodons() : CCSignal() {}
  virtual ~CCSignalWithCodons() {}

  virtual void SetSequence(const DNAVector & b, int size);

  virtual int GetFullSize() const {return GetA().size();}
  
  virtual void SetSize(int n);

  virtual int GetCount() const {return 4 + 21;}
  virtual const svec<float> & Get(int i) const {
    if (i < 4)
      return CCSignal::Get(i);
    return m_aa[i-4];
  }

 private:
  vecFloat m_aa[21];
};






class CrossCorrelation
{
 public:
  CrossCorrelation();
  ~CrossCorrelation();
  void CrossCorrelate(svec<float> & out, const CCSignal & one, const CCSignal & two);
  
  void AutoCorrelate(vector<float> &out, vector<float> &in);

  void DoOne(svec<float> & o, const svec<float> & in1, const svec<float> & in2);

 private:
  FFTReal<float> * m_pFFT;
  int m_size;

};



class SeqMatch
{
 public:
  SeqMatch() {
    m_startTarget = -1;
    m_startQuery = -1;
    m_len = 0;
    m_ident = 0.;
  }
 
  SeqMatch(int startTarget, int startQuery, int len, double ident) {
    Set(startTarget, startQuery, len, ident);
  }

  void Set(int startTarget, int startQuery, int len, double ident) {
    m_startTarget = startTarget;
    m_startQuery = startQuery;
    m_len = len;
    m_ident = (float)ident;
  }

  int GetStartTarget() const {return m_startTarget;}
  int GetStartQuery() const  {return m_startQuery;}
  int GetLength() const {return m_len;}
  double GetIdentity() const {return m_ident;}

private:
  int m_startTarget;
  int m_startQuery;
  int m_len;
  float m_ident;
};


class vecSeqMatch
{
 public:
  vecSeqMatch() {
    m_len = 0;
  }

  void clear() {m_len = 0;}

  void push_back(const SeqMatch & m) {
    if (m_len >= m_data.isize())
      m_data.resize(m_len + 4096);
    m_data[m_len] = m;
    m_len++;
  }

  int size() const {return m_len;}
  int isize() const {return m_len;}

  SeqMatch & operator[] (int i) {return m_data[i];}
  const SeqMatch & operator[](int i) const {return m_data[i];}


 private:
  svec<SeqMatch> m_data;
  int m_len;
};


class LookupMatch
{
public:
  LookupMatch();

  int GetScale() const {return m_scale;}

  int GetScore(char a, char b) const {
    int i = m_index[(int)a];
    int j = m_index[(int)b];
    return m_score[16*i + j];
  }

private:
  void Set(char c, int index) {
    m_index[(int)c] = index;
    m_letter[index] = c;
  }

  int m_index[256];
  char m_letter[256];
  int m_scale;
  
  int m_score[256];
};



class SeqAnalyzer
{
 public:
  SeqAnalyzer();
  void MatchUp(vecSeqMatch & out, const DNAVector & query, const DNAVector & target, svec<float> & xc);

  void MatchUp(vecSeqMatch & out, const CCSignal & query, const CCSignal & target, svec<float> & xc);

  void SetTopCutoff(double c) {
    m_topCutoff = c;
    //cout << "Using cutoff " << m_topCutoff << endl;
  }
  void SetMinLen(int l) {m_minLen = l;}
  void SetMinIdent(double m) {m_minIdent = m;}


 private:
  int FindBest(svec<float> & cx);
  int FindTop(svec<int> & top, svec<float> & xc, double topCutoff);
  void DoOne(vecSeqMatch & out, const DNAVector & query, const DNAVector & target, int pos);
  void DoOneFloat(vecSeqMatch & out, const CCSignal & query, const CCSignal & target, int pos);
  void DoOneSlow(vecSeqMatch & out, const DNAVector & query, const DNAVector & target, int pos);
  void DoOneFull(vecSeqMatch & out, const DNAVector & query, const DNAVector & target, int pos);
  bool IsGood(const DNAVector & query, const DNAVector & target, int posQ, int posT, int len);

  double m_minIdent;
  double m_topCutoff;
  int m_minLen;

  svec<double> m_envelope;
  int m_envSize;

  LookupMatch m_lookup;
};


double PrintMatch(const DNAVector & query, const DNAVector & target, const SeqMatch & m, bool bSilent = false, bool bNukes = true);



#endif

