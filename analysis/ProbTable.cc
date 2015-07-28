#include "analysis/ProbTable.h"

ProbTable::ProbTable() {
  m_size = 100000;
  m_cutoff = 0.999;
}
ProbTable::ProbTable(double targetSize, double cutoff) {
  Setup(targetSize, cutoff);
}

void ProbTable::Setup(double targetSize, double cutoff) {
  m_size = targetSize;
  m_cutoff = cutoff;

  // Accurate parameters, but takes time building the table.
  m_table.resize(512);
  int maxLen = 2048;

  // Less accurate, but faster table build.
  //m_table.resize(128);
  //int maxLen = 1024;
  int i, j;

  for (i=1; i<m_table.isize(); i++) {
    svec< double > & t = m_table[i];    
    double ident_expect = (double)i/((double)m_table.isize()-1);
    t.resize(maxLen, 2.);
    for (j=1; j<maxLen; j++) {
      for (double ident = 0.001; ident <=1.; ident += 0.001) {
        double prob = GetMatchProbabilityRaw(j, ident, ident_expect, m_size);
        if (prob > m_cutoff) {
          t[j] = ident;
          break;
        }
      }   
    }
  }
}

bool ProbTable::IsGood(int length, double ident, double ident_expect) {
  int index = ExpectToIndex(ident_expect);
  svec< double > & t = m_table[index];    
  if (length >= t.isize())
    length = t.isize() - 1;
  double cutoff = t[length];
  if (ident >= cutoff) {
    return true;
  }
  return false;
}

double ProbTable::GetMatchProbabilityRaw(int length, double ident, double ident_expect, double targetSize);

// Funtion compatible with the old match
double ProbTable::GetMatchProbability(double & ident,
    const DNAVector & target, 
    const DNAVector & query, 
    int startTarget,
    int startQuery,
    int length);
int ProbTable::ExpectToIndex(double ident_expect) {
  int index = ident_expect * (m_table.isize()-1);
  return index;
}

double ProbTable::CDF(double x, double m, double s)
{

  double r = 0.5 * (1. + erf((x-m)/s/1.414213562));

  return r;
}

double ProbTable::Sigma(double p, int N)
{
  return sqrt(p * (1. - p) * (double)N);
}


double ProbTable::GCAdjustExpect(double gc, int N, double gc_target)
{
  double at_target = 1. - gc_target;

  double r = gc * gc_target;


  r += (double)(N - gc) * at_target;


  return r / (double)N / 2.;

}


double ProbTable::GetMatchProbability(double & ident,
    const DNAVector & target, 
    const DNAVector & query, 
    int startTarget,
    int startQuery,
    int length)
{

  int i;
  double gcCountTarget = 0;
  double gcCountQuery = 0;
  double matches = 0;

  int trimLeft = -1;
  int lastMatch = 0;

  for (i=0; i<length; i++) {
    if (target[i+startTarget] == query[i+startQuery]) {
      if (trimLeft == -1)
        trimLeft = i;
      lastMatch = i;
    }

    matches += DNA_Equal(target[i+startTarget], query[i+startQuery]);

    gcCountTarget += DNA_C(target[i+startTarget]) + DNA_G(target[i+startTarget]); 
    gcCountQuery += DNA_C(query[i+startQuery]) +  DNA_G(query[i+startQuery]);

  }

  ident = matches/(double)length;

  double gc = gcCountQuery;

  double p_match = GCAdjustExpect(gc, length, gcCountTarget/(double)length);
  bool b = IsGood(length, ident, p_match);
  if (b) {
    return m_cutoff;
  } else {
    return 0.;
  }
}


double ProbTable::GetMatchProbabilityRaw(int length, double ident, double ident_expect, double targetSize)
{

  int i;


  double p_match = ident_expect;



  double s = Sigma(p_match, length);

  double m = p_match * (double)length;

  double x = (double)length * ident;

  double cdf = CDF(m, x, s);
  double expect = cdf * targetSize; // Twice the genome size if we do RC!!
  double p_val = exp(-expect);


  return p_val;
}
