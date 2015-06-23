#ifndef FORCE_DEBUG
#define NDEBUG
#endif



 
#include "analysis/CrossCorr.h"
#include "analysis/DNAVector.h"
#include "analysis/CodonTranslate.h"
#include "extern/RealFFT/FFTReal.hpp"
#include "analysis/Blosum.h"
#include <math.h>






void ComplexMult(float & r_o, float & i_o, float r1, float i1, float r2, float i2) 
{
  r_o = r1 * r2 - i1 * i2;
  i_o = i1 * r2 + i2 * r1;

  //  cout << "r1=" << r1 << " i1=" << i1 << " r2=" << r2 << " i2=" << i2 << " r_o=" << r_o << " i_o=" << i_o << endl;
}

double Ent(double p) 
{
  if (p < 0.001)
    return 0;
  return p * log(p) / 0.69314718056;
}

void CCSignal::ComputeEntropy(const DNAVector & in)
{
  int i, j;

  if (in.isize() < 1024) {
    //cout << "Skip entropy." << endl
    for (i=0; i<m_entropy.isize(); i++)
      m_entropy[i] = 1.;
    return;
  }

  int win = m_entropy.isize() / 512; // 32 bp windows
  for (i=0; i<(int)in.size(); i+= win) {
    double a = 0;
    double c = 0;
    double g = 0;
    double t = 0;
    int k = 0;

    for (j=i; j<i+win; j++) {
      if (j >= (int)in.size())
	break;
      k++;

      /*
      switch(in[j]) {
      case 'A':
	a += 1.;
	break;
      case 'C':
	c += 1.;
	break;
      case 'G':
	g += 1.;
	break;
      case 'T':
	t += 1.;
	break;
	}*/
      a += DNA_A(in[j]);
      c += DNA_C(in[j]);
      g += DNA_G(in[j]);
      t += DNA_T(in[j]);
    }

    a /= (double)k;
    c /= (double)k;
    g /= (double)k;
    t /= (double)k;
    double s = Ent(a) + Ent(c) + Ent(g) + Ent(t);
    for (j=i; j<i+k; j++) {
      m_entropy[j] = -s;
      if (m_entropy[j] < 0.)
	m_entropy[j] = 0.;
    }
    //cout << "i=" << i << " entropy=" << m_entropy[i] << endl;
  }

}



void CCSignal::SeqToPCM(svec<float> & out, const DNAVector & in, char nuke) 
{
  //========================================================
  int i, j;

//=========================================================
 
  int s = m_points;


  // First, compute the average.
  double sum = 0;
  for (i=0; i<(int)in.size(); i++) {
    sum += DNA_N(in[i], nuke);
    //char b = in[i];
    //if (b == nuke)
    //sum++;
  }

  double off = sum / (double)in.size();
  
  int k = 0;
  for (i=0; i<(int)in.size(); i++) {
    char b = in[i];

    //cout << "Setting  nuke " << nuke << " t=" << b << " " << DNA_N(b, nuke) << " ent=" << m_entropy[k] << endl;
    out[k] = m_entropy[k] * (DNA_N(b, nuke) - off);
    //cout << " val=" << out[k] << endl;

    //if (b == nuke)     
    //  out[k] = (1. - off) * m_entropy[k];
    //else
    //  out[k] = - off * m_entropy[k];
    k++;
 
  }
 
}

void Smooth(svec<float> & out, const svec<float> & in)
{
  int i;
  out.clear();
  out.resize(in.isize(), 0);
  for (i=1; i<in.isize()-1; i++) {
    float d = in[i] + (in[i-1] + in[i+1])/2;
    out[i] = d / 2;
  }
}




CCSignal::CCSignal()
{
  m_points = 1;
}

void CCSignal::SetSize(int n)
{ 
  int s = m_points;

  // 
  int len = n * s * 2;

  //cout << "Real size: " << len << endl;
  int i;

  int b = 1;
  while (b < len) {
    b *= 2;
  }
  len = b;
  //cout << "Assigned size=" << len << endl;

  m_A.resize(len, 0.);
  m_C.resize(len, 0.);
  m_G.resize(len, 0.);
  m_T.resize(len, 0.);
  m_entropy.resize(len, 0.);
}

void CCSignal::SetSequence(const DNAVector & b, int size)
{

  if (size != m_A.isize()) {
    m_A.clear();
    m_C.clear();
    m_G.clear();
    m_T.clear();
    
    if (size == 0)
      SetSize((int)b.size());
    else {
      m_A.resize(size, 0.);
      m_C.resize(size, 0.);
      m_G.resize(size, 0.);
      m_T.resize(size, 0.);
      m_entropy.resize(size, 0.);
    }
  }

  ComputeEntropy(b);

  SeqToPCM(m_A, b, 'A');
  SeqToPCM(m_C, b, 'C');
  SeqToPCM(m_G, b, 'G');
  SeqToPCM(m_T, b, 'T');

}


void CCSignal::Smooth()
{
  svec<float> tmp;
  ::Smooth(tmp, m_A);
  ::Smooth(m_A, tmp);
  ::Smooth(tmp, m_C);
  ::Smooth(m_C, tmp);
  ::Smooth(tmp, m_G);
  ::Smooth(m_G, tmp);
  ::Smooth(tmp, m_T);
  ::Smooth(m_T, tmp);  
}


void CCSignalProtein::SetSequence(const DNAVector & b, int size)
{
  int i, j;
  if (size != m_aa[0].isize()) {
    for (i=0; i<21; i++) {
      m_aa[i].clear();
      m_aa[i].resize(size, 0.);
    }
  }
  


  for (i=0; i<21; i++) {
    double sum = 0;
    char ll = AminoAcidLetter((char)i);
    vecFloat & aa = m_aa[i];
    for (j=0; j<(int)b.size(); j++) { 
      int sim = Similarity(ll, b[j]);
      double add = 0.;
      if (sim > 0)
	add = ((double)sim)/30; // Hard coded!!!
      if (sim == -100) {
	//cout << "WARNING: Amino acid not found in Blosum table: " << b[j] << endl;
	add = 0.;
      }
      if (b[j] == ll) {
	aa[j] = 1.;
 	sum += 1.;
      } else {
	//aa[j] = add;
	//sum += add;
      }
    }
    
    double off = sum / (double)b.size();
  
    for (j=0; j<(int)b.size(); j++) { 
      aa[j] -= off;
    }
    for (; j<aa.isize(); j++)
      aa[j] = 0.;
    
    double check = 0;
    for (j=0; j<aa.isize(); j++)
      check += aa[j];
    //cout << "Check: " << check << endl;
  }
  

}

void CCSignalProtein::SetSize(int n)
{
  CCSignal::SetSize(n);
  for (int i=0; i<21; i++) {
    m_aa[i].clear();
    m_aa[i].resize(GetFullSize(), 0.);
  }
}

//==================================================================
void CCSignalWithCodons::SetSequence(const DNAVector & b, int size)
{
  int i;
  CCSignal::SetSequence(b, size);
  if (size != 0) {
    for (i=0; i<21; i++) {
      m_aa[i].clear();
      m_aa[i].resize(size, 0.);
    }
  }
 
  static CodonTranslator trans;
  
  cout << "Start w/ base " << b[0] << b[1] << b[2] << endl;

  for (i=0; i<(int)b.size()-2; i+= 3) {
    char b1 = b[i];
    char b2 = b[i+1];
    char b3 = b[i+2];
    int idx = trans.GetCodonIndex(b1, b2, b3);
    //cout << b1 << b2 << b3 << "  idx=" << idx << endl;

    const CodonMatrix & m = trans.GetMatrix(idx);
    
    for (int j=0; j<m.Size(); j++) {
      svec<float> & f = m_aa[j];
      f[i] = m.Get(j);
      f[i+1] = m.Get(j);
      f[i+2] = m.Get(j);
    }

    //svec<float> & f2 = m_aa[idx];
    //f2[i] = 1.;

  }

  int j;

  
  for (i=0; i<21; i++) {
    svec<float> & f = m_aa[i];
    double sum = 0.;
    for (j=0; j<(int)b.size(); j++) 
      sum += f[j];
    
    //cout << "Protein # " << i << endl;
    //cout << "sum=" << sum << endl;
    sum /= (double)b.size();
    for (j=0; j<(int)b.size(); j++) {
      f[j] -= sum;
      //cout << "j=" << j << " val=" << f[j] << endl;
    }
  }
}


void CCSignalWithCodons::SetSize(int n)
{
  CCSignal::SetSize(n);
  for (int i=0; i<21; i++) {
    m_aa[i].clear();
    m_aa[i].resize(GetFullSize(), 0.);
  }
}

//==================================================================
CrossCorrelation::CrossCorrelation()
{
  m_pFFT = NULL;
  m_size = 0;
}

CrossCorrelation::~CrossCorrelation()
{
  if (m_pFFT != NULL)
    delete m_pFFT;
}


void CrossCorrelation::AutoCorrelate(vector<float> &out, vector<float> &in)
{
  out.clear();

  int N=FFTReal_get_next_pow2((long) in.size());

  in.resize(N,0);
  out.resize(N,0);
  
  vector<float> X(N,0),tmp(N,0);

  m_pFFT = new FFTReal<float> (N);

  m_pFFT->do_fft(&X[0],&in[0]);
  
  for (int i=0; i<=N/2; ++i )
    tmp[i] = (X[i]+X[i+N/2])*(X[i]-X[i+N/2]);

  m_pFFT->do_ifft(&tmp[0],&out[0]);
  m_pFFT->rescale(&out[0]);
}


void CrossCorrelation::CrossCorrelate(svec<float> & out, const CCSignal & one, const CCSignal & two)
{
  out.clear();
  out.resize(one.GetFullSize(), 0.);

  svec<float> tmp;
  int i, j;

 
  for (j=0; j<one.GetCount(); j++) {
    DoOne(tmp, one.Get(j), two.Get(j));

    for (i=0; i<out.isize(); i++) {
      //if (one.GetCount() <= 4 || j > 4)
      out[i] += tmp[i];
      tmp[i] = 0.;
    }
  }


  /*
  DoOne(tmp, one.GetA(), two.GetA());

  for (i=0; i<out.isize(); i++) {
    out[i] += tmp[i];
    //out[i] = tmp[i];
    tmp[i] = 0.;
  }
  
  DoOne(tmp, one.GetC(), two.GetC());
  
  for (i=0; i<out.isize(); i++) {
    out[i] += tmp[i];
    //out[i] *= tmp[i];
    tmp[i] = 0.;
  }
  
  DoOne(tmp, one.GetG(), two.GetG());
  
  for (i=0; i<out.isize(); i++) {
    out[i] += tmp[i];
    //out[i] *= tmp[i];
    tmp[i] = 0.;
  }
  DoOne(tmp, one.GetT(), two.GetT());
  

  for (i=0; i<out.isize(); i++) {
    out[i] += tmp[i];
    //out[i] *= tmp[i];
    //cout << i << " -> " << out[i] << endl;
    tmp[i] = 0.;
  }
  */

}



void CrossCorrelation::DoOne(svec<float> & o, const svec<float> & in1, const svec<float> & in2)
{
  int i;
  o.resize(in1.size(), 0);

  if (m_pFFT == NULL || m_size != in1.isize()) {
    m_size = in1.isize();
    if (m_pFFT != NULL) {
      cout << "WARNING: re-instantiating FFT object!" << endl;
      delete m_pFFT;
    }
    //cout << "Initializing FFT." << endl;
    m_pFFT = new FFTReal<float>(m_size);
    //cout << "done." << endl;
  }

  svec<float> tmp1;
  tmp1.resize(in1.size(), 0.);
  svec<float> tmp2;
  tmp2.resize(in2.size(), 0.);

  float * p1 = &tmp1[0];
  float * p2 = &tmp2[0];



  m_pFFT->do_fft(p1, &in1[0]);


  m_pFFT->do_fft(p2, &in2[0]);


  int N = tmp1.isize() / 2;
  //cout << "N=" << N << endl;

  tmp1[0] *= tmp2[0];
  //tmp1[N] *= tmp2[N];
  for (i=1; i<N-1; i++) {
    float fr, fi;

    //ComplexMult(fr, fi, tmp1[i], tmp1[i+N+1], tmp2[i], -tmp2[i+N+1]);
    //tmp1[i] = fr;
    //tmp1[i+N+1] = -fi;

    ComplexMult(fr, fi, tmp1[i], tmp1[i+N], tmp2[i], -tmp2[i+N]);
    tmp1[i] = fr;
    tmp1[i+N] = -fi;
  }
  m_pFFT->do_ifft(p1, p2);
  m_pFFT->rescale(p2);

  //for (i=0; i<2*N; i++)
  //o[i] = tmp2[i];

  
  for (i=0; i<N; i++) {
    o[i] = tmp2[i+N];
  }
  for (i=N; i<2*N; i++) {
    o[i] = tmp2[i-N];
  }
 
}




//============================================================
LookupMatch::LookupMatch() { 
  int i, j;
  for (i=0; i<256; i++) {
    m_index[i] = -1;
    m_letter[i] = 0;
  }
  
  Set('A', 0);
  Set('C', 1);
  Set('G', 2);
  Set('T', 3);
  Set('K', 4);
  Set('M', 5);
  Set('R', 6);
  Set('Y', 7);
  Set('S', 8);
  Set('W', 9);
  Set('B', 10);
  Set('V', 11);
  Set('H', 12);
  Set('D', 13);
  Set('N', 14);
  Set('X', 15);
  
  m_scale = 100;
  //m_score.resize(256);
  
  for (i=0; i<16; i++) { 
    for (j=0; j<16; j++) { 
      double d = DNA_EqualAmb(m_letter[i], m_letter[j]);
      int v = (int)(d * (double)m_scale + 0.5);
      m_score[16 * i + j] = v;
    }
  }
}


SeqAnalyzer::SeqAnalyzer()
{
  m_minIdent = 0.42;
  m_minLen = 45;
  m_topCutoff = 1.8;
  m_envSize = 256;
}


void SeqAnalyzer::MatchUp(vecSeqMatch & out, const CCSignal & query, const CCSignal & target, svec<float> & xc)
{
  int i;

  svec<int> all;
  all.reserve(400);

  FindTop(all, xc, m_topCutoff);

  for (i=0; i<all.isize(); i++) {
    int pos = all[i];
    pos -= xc.isize() / 2;
    DoOneFloat(out, query, target, pos);
  }
}

 
void SeqAnalyzer::MatchUp(vecSeqMatch & out, const DNAVector & query, const DNAVector & target, svec<float> & xc)
{
  int i;
  
  svec<int> all;
  all.reserve(400);

  FindTop(all, xc, m_topCutoff);

  //all.resize(xc.isize());
  //for (i=0; i<xc.isize(); i++) {
  //  all[i] = i;
  //}

  //cout << "Evaluating " << all.isize() << " positions." << endl;

  for (i=0; i<all.isize(); i++) {
    int pos = all[i];
    pos -= xc.isize() / 2;
    DoOne(out, query, target, pos);
    //DoOneSlow(out, query, target, pos);
  }
}



#define OPTIM_DOONE
#ifdef OPTIM_DOONE
void SeqAnalyzer::DoOneFloat(vecSeqMatch & out, const CCSignal & query, const CCSignal & target, int pos)
{
  int i, j;
  int shift = pos;

  double match = 0;
  
  double minMatches = (int)((double)m_minLen * m_minIdent);

  int lastStart = -1;

  int n = 0;

  int minLen = m_minLen;
  if (target.GetFullSize() + 4 < minLen)
    minLen = target.GetFullSize()-4;

  for (i=0; i<(int)target.GetFullSize(); i++) {
  
    int j = i + shift;

    if (j < 0)
      continue;
    if (j >= (int)query.GetFullSize() || i + 1 >= (int)target.GetFullSize()) {
      if (lastStart != -1) {
	out.push_back(SeqMatch(lastStart, lastStart + shift, i-lastStart, 0.));
      }
      break;
    } 

    //match += DNA_EqualFast(target[i], query[j]);
    match += CCScore(target, query, i, j);
    //cout << " matchplus " << i << "\t" << match << endl;
    if (n > m_minLen) {
   
      //match -= DNA_EqualFast(target[i-m_minLen-1], query[j-m_minLen-1]);
      match -= CCScore(target, query, i-m_minLen-1, j-m_minLen-1);

      //cout << " match " << i << "\t" << match << endl;

      if (match > minMatches) {
	if (lastStart == -1)
	  lastStart = i - m_minLen;
      } else {
	if (lastStart != -1) {
	  out.push_back(SeqMatch(lastStart, lastStart + shift, i-lastStart, 0.));
	}
	lastStart = -1;	
      }
    }

    n++;

  }
}

void SeqAnalyzer::DoOne(vecSeqMatch & out, const DNAVector & query, const DNAVector & target, int pos)
{
  int i, j;
  int shift = pos;


  int match = 0;
  
  int minMatches = (int)((double)m_minLen * m_minIdent * (double)m_lookup.GetScale());

  int lastStart = -1;

  int n = 0;

  int minLen = m_minLen;
  if (target.isize() + 4 < minLen)
    minLen = target.isize()-4;

  for (i=0; i<(int)target.size(); i++) {
  
    int j = i + shift;

    if (j < 0)
      continue;
    if (j >= (int)query.size() || i + 1 >= (int)target.size()) {
      if (lastStart != -1) {
	out.push_back(SeqMatch(lastStart, lastStart + shift, i-lastStart, 0.));
      }
      break;
    } 

    //match += DNA_EqualFast(target[i], query[j]);
    match += m_lookup.GetScore(target[i], query[j]);
    if (n > m_minLen) {
   
      //match -= DNA_EqualFast(target[i-m_minLen-1], query[j-m_minLen-1]);
      match -= m_lookup.GetScore(target[i-m_minLen-1], query[j-m_minLen-1]);


      if (match > minMatches) {
	if (lastStart == -1)
	  lastStart = i - m_minLen;
      } else {
	if (lastStart != -1) {
	  out.push_back(SeqMatch(lastStart, lastStart + shift, i-lastStart, 0.));
	}
	lastStart = -1;	
      }
    }

    n++;

  }
  
  //cout << "Done here!" << endl;
}

#else //=================================================================

void SeqAnalyzer::DoOne(vecSeqMatch & out, const DNAVector & query, const DNAVector & target, int pos)
{
  int i, j;
  int shift = pos;


  double match = 0.;
  
  double minMatches = (int)((double)m_minLen * m_minIdent);

  int lastStart = -1;

  int n = 0;


  for (i=0; i<(int)target.size(); i++) {
  
    int j = i + shift;

    if (j < 0)
      continue;
    if (j >= (int)query.size() || i + 1 >= (int)target.size()) {
      if (lastStart != -1) {
	out.push_back(SeqMatch(lastStart, lastStart + shift, i-lastStart, 0.));
      }
      break;
    } 

    match += DNA_EqualFast(target[i], query[j]);

    if (n > m_minLen) {
   
      match -= DNA_EqualFast(target[i-m_minLen-1], query[j-m_minLen-1]);


      if (match > minMatches) {
	if (lastStart == -1)
	  lastStart = i - m_minLen;
      } else {
	if (lastStart != -1) {
	  out.push_back(SeqMatch(lastStart, lastStart + shift, i-lastStart, 0.));
	}
	lastStart = -1;	
      }
    }

    n++;

  }
  
  //cout << "Done here!" << endl;
}
#endif
//==========================================================================


void SeqAnalyzer::DoOneFull(vecSeqMatch & out, const DNAVector & query, const DNAVector & target, int pos)
{
  int tLen = (int)target.size();
  int qLen = (int)query.size();
  
  //if (tLen != qLen) {
  //cout << "SKIP!" << endl;
  //return;
  //}

  if (pos >= 0) {
    int lapLen = qLen - pos;
    if (tLen < lapLen)
      lapLen = tLen;
    cout << "PUSH + " << pos << " " << lapLen << endl;
    out.push_back(SeqMatch(0, pos, lapLen, 0.));
    //out.push_back(SeqMatch(pos, 0, lapLen, 0.));
  } else {
    int lapLen = tLen + pos;
    if (qLen < lapLen)
      lapLen = qLen;    
    cout << "PUSH - " << -pos << " " << lapLen << endl;
    out.push_back(SeqMatch(-pos, 0, lapLen, 0.));
    //out.push_back(SeqMatch(0, -pos, lapLen, 0.));
  }
  
}

void SeqAnalyzer::DoOneSlow(vecSeqMatch & out, const DNAVector & query, const DNAVector & target, int pos)
{
  int i, j;
  int shift = pos;

  //cout << "Shifting sequences by " << shift << " bases." << endl;

  double match = 0;
  
  //svec<int> matches;
  //matches.resize(target.size(), 0);
  double minMatches = ((double)m_minLen * m_minIdent);

  int lastStart = -1;
  for (i=0; i<(int)target.size()-m_minLen; i++) {
    int n = 0;
    match = 0.;
    //cout << "top=" << i + m_minLen << endl;
    for (j=i; j<i+m_minLen; j++) {
      int qi = j + shift;

      if (qi < 0 || qi >= (int)query.size())
	continue;
      n++;

      
      //if (target[j] == query[qi]) {
      //	match++;
      //}
      match += DNA_EqualEmph(target[j], query[qi]);
    }
    if (n == 0)
      continue;
    if (match < minMatches || n < m_minLen || j+1 == (int)target.size()) {
      if (lastStart != -1) {
	//cout << "Adding match, ls=" << lastStart << " length=" << j-lastStart-1 << endl;
	//cout << "match=" << match << " n=" << n << endl;
	if (IsGood(query, target, lastStart+shift, lastStart, j-lastStart-1))
	  out.push_back(SeqMatch(lastStart, lastStart + shift, j-lastStart-1, 0.));
      }
      lastStart = -1;	
    } else {
      if (n == m_minLen && lastStart == -1)
	lastStart = j - m_minLen;
    }
  }
  //cout << "Done here!" << endl;
}

int SeqAnalyzer::FindBest(svec<float> & xc)
{
  int i;
  double max = 0.;
  int best = -1;
  for (i=0; i<xc.isize(); i++) {
    //cout << "i=" << i << " -> " << xc[i] << endl;
    if (xc[i] > max) {
      max = xc[i];
      best = i;
    }
  }
  //cout << "Best: " << best << "  " << max << endl;

  return best;
}

int SeqAnalyzer::FindTop(svec<int> & top, svec<float> & xc, double topCutoff)
{
  //First, let's find the best score
  int i;
  double max = 0.;
  int best = -1;
  double overTop = 1.;

  int nPoints = xc.isize() / m_envSize;
  if (nPoints == 0)
    nPoints = 1;
  //cout << "nPoints: " << nPoints << endl;
  if (m_envelope.isize() != nPoints) {
    m_envelope.clear();
    m_envelope.resize(nPoints, 0.) ;
  } else {
    for (i=0; i<m_envelope.isize(); i++)
      m_envelope[i] = 0.;
  } 


  double avg = 0.;

  for (i=0; i<xc.isize(); i++) {
    double d = xc[i] * xc[i];
    avg += xc[i];
    //cout << "i=" << i << " index=" << i/nPoints << endl;
    if (nPoints > 8)
      m_envelope[i/m_envSize] += d;
    if (xc[i] > max) {
      max = xc[i];
      best = i;
    }
  }

  // MGG: Possible speed-up: exit if the average is below an absolute  threshold
  avg /= (double)xc.isize();
  //cout << "Average=" << avg << endl;
  //if (avg < 0.005) {
  //  top.push_back(best);
  //  return best;
  // }


  for (i=0; i<m_envelope.isize(); i++) {
    m_envelope[i] = sqrt(m_envelope[i] / (double)m_envSize);
    //cout << "i=" << i << " env=" << m_envelope[i] << endl;
  }
  

  max *= topCutoff;
  // Get everything in range
  for (i=0; i<xc.isize(); i++) {
    //cout << "i=" << i << " xc=" << xc[i] << " dev=" << m_envelope[i/m_envSize];
    if (xc[i] > m_envelope[i/m_envSize] * topCutoff + overTop) {
      //cout << "Accepting " << xc[i] << " " << i << " over " <<  m_envelope[i/m_envSize] * topCutoff + overTop << endl;
      top.push_back(i);
      //cout << " ****";
    }
    //cout << endl;
    //if (xc[i] > max) {
    //top.push_back(i);    
    //}
  }

  return best;
}

bool SeqAnalyzer::IsGood(const DNAVector & query, const DNAVector & target, int posQ, int posT, int len)
{
  return true;


  int i;
  int match = 0;
  for (i=0; i<len; i++) {
    if (query[i+posQ] == target[i+posT])
      match++;
  }

  double expect = (double)len / 3.;
  double sigma = sqrt(expect);
  double minAccept = expect + 2 * sigma;
  if ((double)match > minAccept)
    return true;
  else
    return false;  
}

double PrintMatch(const DNAVector & query, const DNAVector & target, const SeqMatch & m, bool bSilent, bool bNuc)
{
  int i, j;

  svec<string> q, t, a;
  int l = 80;
 
  string qq, tt, mm;
  int k = 0;

  //cout << m.GetStartTarget() << "  " <<  m.GetLength() << endl;
  double matches = 0;
  for (i=m.GetStartTarget(); i<m.GetStartTarget() + m.GetLength(); i++) {
    if ((k > 0 && k % l == 0) || i+1 == m.GetStartTarget() + m.GetLength()) {
      q.push_back(qq);
      t.push_back(tt);
      a.push_back(mm);
      qq = "";
      tt = "";
      mm = "";
    }
    j = m.GetStartQuery() + i - m.GetStartTarget();
    //cout << "i=" << i << " j=" << j << endl;
    qq += query[j];
    tt += target[i];

    double dist = 0;
    if (bNuc) {
      dist = DNA_Equal(query[j], target[i]);
    } else {
      if (query[j] == target[i])
	dist = 1.;
    }
    matches += dist;
    
    string matchInd = " ";
    if (dist > 0.3) {
      matchInd = ".";
    }
    if (dist > 0.49) {
      matchInd = "~";
    }
    if (dist > 0.9) {
      matchInd = "|";
    }

    mm += matchInd;

    //if (query[j] == target[i]) {
    //  mm += "|";
    //} else {
    //  mm += " ";
    //}
    k++;
  }

  if (!bSilent) {
    cout << "Identity = " << 100 * matches/(double)m.GetLength() << " % ";
    cout << "  Alignment length = " << m.GetLength() << endl;
    cout << "Query " << m.GetStartQuery() << " - " << m.GetStartQuery() + m.GetLength();
    cout << "   Target " << m.GetStartTarget() << " - " << m.GetStartTarget() + m.GetLength() << endl;
    cout <<  "--------------------------------------------------------------------------------" << endl;      
    for (i=0; i<q.isize(); i++) {
      cout << q[i] << endl;
      cout << a[i] << endl;
      cout << t[i] << endl;
      cout << endl << "--------------------------------------------------------------------------------" << endl << endl;      
    }
  }
  return (double)matches/(double)m.GetLength();
}

