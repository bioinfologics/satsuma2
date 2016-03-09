#include "ProbTable.h"

ProbTable::ProbTable() {
  m_size = 100000;
  m_cutoff = 0.999;
}
ProbTable::ProbTable(double targetSize, double cutoff) {
  Setup(targetSize, cutoff);
  for (unsigned char i=0;i<128;i++){
    for (unsigned char j=0;j<128;j++) dna_id_table[i][j]=DNA_Equal(i,j);
    dna_gc_table[i]= DNA_C(i) + DNA_G(i); 
  }
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

  for (i=1; i<m_table.size(); i++) {
    std::vector< double > & t = m_table[i];    
    double ident_expect = (double)i/((double)m_table.size()-1);
    //cout<<"Probtable ident_expect="<<ident_expect<<endl;
    t.resize(maxLen, 2.);
    double ident;
    double min,max,mid,prob;
    for (j=1; j<maxLen; j++) { //XXX you can do a freaking binary search for the limit!
      for (min=0,max=1;max-min>0.00000001;){
        mid=(max+min)/2.0;
        if (GetMatchProbabilityRaw(j, mid, ident_expect, m_size))max=mid;
        else min=mid;
      }
      t[j]=(max+min)/2.0;
    /*
      //binary search with 10 digits
      for (ident = 0.001; ident <=1.; ident += 0.001) {
        double prob = GetMatchProbabilityRaw(j, ident, ident_expect, m_size);
        //cout<<"GetMatchProbabilityRaw(j="<<j<<", ident="<<ident<<", ident_expect="<<ident_expect<<", m_size="<<m_size<<") = "<<prob<<endl;
        if (prob > m_cutoff) {
          t[j] = ident;
          //cout<<"ProbTable["<<i<<"]["<<j<<"] = "<<ident<<endl;
          break;
        }
      }   
      //if (ident >1.) cout<<"ProbTable["<<i<<"]["<<j<<"] = "<<t[j]<<" (defaulting)"<<endl;*/
    }
  }
}

bool ProbTable::IsGood(int length, double ident, double ident_expect) {
  int index = ExpectToIndex(ident_expect);
  std::vector< double > & t = m_table[index];    
  if (length >= t.size())
    length = t.size() - 1;
  double cutoff = t[length];
  //std::cout<<"Cutoff for length "<< length<<" is "<<cutoff<<std::endl;
  if (ident >= cutoff) {
    return true;
  }
  return false;
}

int ProbTable::ExpectToIndex(double ident_expect) {
  int index = ident_expect * (m_table.size()-1);
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

  //XXX: optimize this, make sure it can be vectorised.
  for (i=0; i<length; i++) {
    matches += dna_id_table[target[i+startTarget]][query[i+startQuery]];
    gcCountTarget += dna_gc_table[target[i+startTarget]]; 
    gcCountQuery += dna_gc_table[query[i+startQuery]];
    //XXX: change this to calculate GC /(gcat) (copes with Ns a lot better).

  }
  //change it to a sum of lookups id_table[s1[start1+i]][s2[start1+i]]

  ident = matches/(double)length;

  double gc = gcCountQuery;
  //another lookup table?
  double p_match = GCAdjustExpect(gc, length, gcCountTarget/(double)length);
  //std::cout<<"probtable: p_match = "<<p_match<<"  length = "<<length<<"  ident = "<<ident<<endl;
  bool b = IsGood(length, ident, p_match);
  if (b) {
    return m_cutoff;//Do we want to recalculate probability here? it may not be too expensive after all
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
